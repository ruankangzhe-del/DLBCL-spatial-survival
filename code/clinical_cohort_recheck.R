args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1],winslash='/')
# Use the explicitly installed R runtime/library; versions locked in renv.lock.
out <- normalizePath(args[2],winslash='/',mustWork=FALSE)
for (x in c('outputs','audit','logs')) dir.create(file.path(out,x),recursive=TRUE,showWarnings=FALSE)
library(survival)
library(metafor)
options(digits=17)
cfg <- read.delim(file.path(root,'config/clinical_cohorts.tsv'),stringsAsFactors=FALSE,check.names=FALSE)
read_tsv <- function(path) read.delim(path,stringsAsFactors=FALSE,check.names=FALSE,na.strings=c('','NA'))
write_tsv <- function(x,name,folder='outputs') write.table(x,file.path(out,folder,name),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
effects <- list(); diagnostics <- list(); feasibility <- list(); interactions <- list(); slopes <- list(); checks <- list(); fits <- list(); nested <- list()
safe_fit <- function(form,data) {
  warns <- character()
  fit <- tryCatch(withCallingHandlers(coxph(form,data=data,ties='efron',singular.ok=FALSE,x=TRUE,y=TRUE,model=TRUE,control=coxph.control(eps=1e-11,iter.max=1000)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
  if(inherits(fit,'error')) return(list(ok=FALSE,reason=conditionMessage(fit)))
  ok <- length(warns)==0 && all(is.finite(coef(fit))) && all(is.finite(vcov(fit))) && all(diag(vcov(fit))>0) && max(fit$iter)<1000
  list(ok=ok,fit=fit,reason=if(ok) 'PASS' else paste(c(warns,'nonfinite/covariance/convergence screen failed'),collapse=' | '))
}
add_ph <- function(fit,cohort,label,n,events) {
  z <- tryCatch(cox.zph(fit,transform='rank'),error=function(e)e)
  if(inherits(z,'error')) {
    diagnostics[[length(diagnostics)+1]] <<- data.frame(cohort=cohort,model=label,n=n,events=events,term='GLOBAL',chisq=NA,df=NA,p=NA,transform='rank',status='NOT_ESTIMABLE',reason=conditionMessage(z));return()
  }
  for(term in rownames(z$table)) diagnostics[[length(diagnostics)+1]] <<- data.frame(cohort=cohort,model=label,n=n,events=events,term=term,chisq=z$table[term,'chisq'],df=z$table[term,'df'],p=z$table[term,'p'],transform='rank',status='ESTIMATED',reason='NA')
}
ipi_levels <- c('Low','Low-Intermediate','High-Intermediate','High')
for(i in seq_len(nrow(cfg))) {
  c <- cfg[i,]; s <- read_tsv(file.path(root,c$score_source)); clin <- read_tsv(file.path(root,c$clinical_source))
  if('model' %in% names(s)) s <- s[s$model=='elastic_net',]
  if('cohort' %in% names(s)) s <- s[s$cohort==c$cohort,]
  stopifnot(!anyDuplicated(s$patient_id),!anyDuplicated(clin$patient_id),setequal(s$patient_id,clin$patient_id))
  # Read clinical fields only from the independently prepared clinical source.
  d <- merge(s[,c('patient_id',c$score_field,c$raw_score_field)],clin,by='patient_id',all=TRUE,sort=TRUE)
  d$score <- d[[c$score_field]]; d$time <- d[[c$time_field]]; d$event <- d[[c$event_field]]
  stopifnot(all(d$event %in% c(0,1)),all(is.finite(d$score)),all(is.finite(d$time)),all(d$time>0))
  scale_error <- max(abs(d$score-d[[c$raw_score_field]]/c$original_development_sd))
  stopifnot(scale_error<1e-12,nrow(d)==c$expected_patients,sum(d$event)==c$expected_events)
  sf <- safe_fit(Surv(time,event)~score,d); stopifnot(sf$ok); f <- sf$fit
  b <- unname(coef(f)); se <- sqrt(unname(vcov(f)[1,1])); hr <- exp(b)
  deltas <- abs(c(hr,b,se)-c(c$reference_hr,c$reference_log_hr,c$reference_se)); stopifnot(all(deltas<1e-6))
  checks[[length(checks)+1]] <- data.frame(cohort=c$cohort,scale_error=scale_error,hr_delta=deltas[1],log_hr_delta=deltas[2],se_delta=deltas[3],tolerance=1e-6,status='PASS')
  effects[[i]] <- data.frame(cohort=c$cohort,role=c$role,include_main_meta=c$include_main_meta,n=nrow(d),events=sum(d$event),log_hr=b,se=se,hr=hr,lower95=exp(b-qnorm(.975)*se),upper95=exp(b+qnorm(.975)*se),p=2*pnorm(-abs(b/se)),score_unit='1 original development SD',original_development_sd=c$original_development_sd,ties='efron')
  fits[[paste0(c$cohort,'_score_only')]] <- f
  add_ph(f,c$cohort,'score_only',nrow(d),sum(d$event))
  if(!c$include_main_meta) next
  if(c$cohort=='GSE125966') {
    stopifnot(all(na.omit(d$ipi) %in% ipi_levels),all(na.omit(d$coo) %in% c('GCB','ABC','UNCLASSIFIED')),all(d$treatment %in% c('Rituximab','GA101')))
    ipi <- d$ipi; coo <- ifelse(d$coo=='UNCLASSIFIED','Unclassified',d$coo); tx_levels <- c('Rituximab','GA101')
  } else {
    stopifnot(all(d$ipi %in% 0:5),all(d$coo %in% c('GCB','ABC','UNC')),all(d$treatment %in% c('R-CHOP','RB-CHOP')))
    ipi <- ipi_levels[findInterval(d$ipi,c(0,2,3,4))]; coo <- ifelse(d$coo=='UNC','Unclassified',d$coo); tx_levels <- c('R-CHOP','RB-CHOP')
  }
  d$ipi4 <- factor(ipi,levels=ipi_levels); d$coo <- factor(coo,levels=c('GCB','ABC','Unclassified')); d$treatment <- factor(d$treatment,levels=tx_levels)
  d <- d[complete.cases(d[,c('score','time','event','ipi4','coo','treatment')]),]
  stopifnot(nrow(d)==c$clinical_n,sum(d$event)==c$clinical_events)
  base_form <- Surv(time,event)~score+ipi4+coo+treatment
  base <- safe_fit(base_form,d)
  if(base$ok) {add_ph(base$fit,c$cohort,'base_adjusted',nrow(d),sum(d$event));fits[[paste0(c$cohort,'_base_adjusted')]]<-base$fit}
  for(v in c('ipi4','coo','treatment')) {
    family <- if(v=='ipi4') 'IPI' else if(v=='coo') 'COO' else 'treatment'
    form <- as.formula(paste('Surv(time,event) ~ score + ipi4 + coo + treatment + score:',v))
    X <- model.matrix(form,d); X <- X[,colnames(X)!='(Intercept)',drop=FALSE]
    ncoef <- ncol(X); rank <- qr(X)$rank
    counts <- table(d[[v]]); deaths <- tapply(d$event,d[[v]],sum)
    reasons <- character()
    if(any(counts<20)) reasons<-c(reasons,'tested_factor_level_patients_below_20')
    if(any(deaths<10)) reasons<-c(reasons,'tested_factor_level_deaths_below_10')
    if(sum(d$event)<10*ncoef) reasons<-c(reasons,paste0('total_deaths_',sum(d$event),'_below_required_',10*ncoef))
    if(rank<ncoef) reasons<-c(reasons,'rank_deficient_design')
    if(!base$ok) reasons<-c(reasons,paste0('base_fit_failed:',base$reason))
    status <- if(length(reasons)) 'NOT_RUN' else 'ESTIMATED'; full <- NULL; lr<-p<-df<-NA_real_
    if(!length(reasons)) {
      full <- safe_fit(form,d)
      if(!full$ok) {status<-'NOT_ESTIMABLE';reasons<-full$reason} else {
        stopifnot(identical(rownames(base$fit$model),rownames(full$fit$model)),base$fit$n==full$fit$n)
        lr <- 2*diff(c(base$fit$loglik[2],full$fit$loglik[2]))
        df <- length(coef(full$fit))-length(coef(base$fit))
        stopifnot(lr>=-1e-8,df==nlevels(d[[v]])-1);lr<-max(0,lr);p<-pchisq(lr,df,lower.tail=FALSE)
        nested[[length(nested)+1]] <- data.frame(cohort=c$cohort,interaction_class=family,base_logLik=base$fit$loglik[2],full_logLik=full$fit$loglik[2],base_coefficients=length(coef(base$fit)),full_coefficients=length(coef(full$fit)),base_n=base$fit$n,full_n=full$fit$n,base_events=base$fit$nevent,full_events=full$fit$nevent,same_rows=identical(rownames(base$fit$model),rownames(full$fit$model)),lr_chisq=lr,df=df,p=p)
        add_ph(full$fit,c$cohort,paste0('score_x_',family),nrow(d),sum(d$event)); fits[[paste0(c$cohort,'_score_x_',family)]]<-full$fit
        for(lev in levels(d[[v]])) {
          contrast <- setNames(rep(0,length(coef(full$fit))),names(coef(full$fit)));contrast['score']<-1
          if(lev!=levels(d[[v]])[1]) {nm<-paste0('score:',v,lev);stopifnot(nm %in% names(contrast));contrast[nm]<-1}
          sl<-sum(contrast*coef(full$fit));slse<-sqrt(drop(t(contrast)%*%vcov(full$fit)%*%contrast))
          slopes[[length(slopes)+1]]<-data.frame(cohort=c$cohort,interaction_class=family,level=lev,n=unname(counts[lev]),events=unname(deaths[lev]),log_hr=sl,se=slse,hr=exp(sl),lower95=exp(sl-qnorm(.975)*slse),upper95=exp(sl+qnorm(.975)*slse),method='linear contrast of common full adjusted model')
        }
      }
    }
    reason <- if(length(reasons)) paste(reasons,collapse=' | ') else 'PASS'
    for(lev in levels(d[[v]])) feasibility[[length(feasibility)+1]] <- data.frame(cohort=c$cohort,interaction_class=family,level=lev,n=unname(counts[lev]),events=unname(deaths[lev]),total_n=nrow(d),total_events=sum(d$event),full_coefficients=ncoef,required_total_events=10*ncoef,design_rank=rank,status=status,reason=reason)
    interactions[[length(interactions)+1]]<-data.frame(cohort=c$cohort,interaction_class=family,n=nrow(d),events=sum(d$event),full_coefficients=ncoef,lr_chisq=lr,df=df,p=p,bh_q=NA_real_,planned_family_n=2,status=status,reason=reason,treatment_contrast=paste(tx_levels[2],'versus',tx_levels[1]),interpretation='exploratory amendment; not original development prespecification')
  }
}
effects <- do.call(rbind,effects); main <- effects[effects$include_main_meta,]; stopifnot(nrow(main)==2)
meta <- rma.uni(yi=main$log_hr,sei=main$se,method='REML',test='adhoc')
pred <- predict(meta)
yy<-main$log_hr; vv<-main$se^2
tau_independent<-max(0,((yy[1]-yy[2])^2-sum(vv))/2)
ww<-1/(vv+tau_independent); mu<-sum(ww*yy)/sum(ww); q<-sum(ww*(yy-mu)^2); qm<-max(1,q); sem<-sqrt(qm/sum(ww));crit<-qt(.975,1)
ind_ci<-mu+c(-1,1)*crit*sem; ind_pi<-mu+c(-1,1)*crit*sqrt(tau_independent+sem^2)
stopifnot(abs(meta$tau2-tau_independent)<1e-10,abs(as.numeric(meta$b)-mu)<1e-10,abs(meta$se-sem)<1e-10,max(abs(c(meta$ci.lb,meta$ci.ub)-ind_ci))<1e-10,max(abs(c(pred$pi.lb,pred$pi.ub)-ind_pi))<1e-10)
meta_out<-data.frame(cohorts=paste(main$cohort,collapse=';'),k=2,method='REML',inference='modified Knapp-Hartung q*=max(1,q)',df=1,log_hr=as.numeric(meta$b),se=meta$se,hr=exp(meta$b[1]),lower95=exp(meta$ci.lb),upper95=exp(meta$ci.ub),p=meta$pval,tau_squared=meta$tau2,I_squared_percent=meta$I2,Q=meta$QE,Q_df=1,Q_p=meta$QEp,q=q,q_modified=qm,prediction_lower95=exp(pred$pi.lb),prediction_upper95=exp(pred$pi.ub),prediction_convention='metafor predict.rma t(df=k-p=1); not Riley df=k-2',limitation='Two cohorts: very uncertain interval; zero estimated heterogeneity does not establish homogeneity')
interactions<-do.call(rbind,interactions)
for(family in unique(interactions$interaction_class)) {idx<-which(interactions$interaction_class==family);interactions$bh_q[idx]<-p.adjust(interactions$p[idx],method='BH',n=2)}
diagnostics <- do.call(rbind,diagnostics)
diagnostics$interpretation <- ifelse(diagnostics$term=='score' & grepl('^score_x_',diagnostics$model),'Score PH at reference level of interacting factor; not all subgroup slopes',ifelse(diagnostics$term=='GLOBAL','Omnibus model PH diagnostic; failure to reject does not prove PH','Term PH diagnostic; failure to reject does not prove PH'))
fit_audit <- do.call(rbind,lapply(names(fits),function(name){f<-fits[[name]];data.frame(model=name,n=f$n,events=f$nevent,free_coefficients=length(coef(f)),design_rank=qr(f$x)$rank,iterations=max(f$iter),log_partial_likelihood=f$loglik[2],finite_coefficients=all(is.finite(coef(f))),finite_covariance=all(is.finite(vcov(f))),minimum_covariance_eigenvalue=min(eigen(vcov(f),symmetric=TRUE,only.values=TRUE)$values),fitting_warning_count=0,ties='efron',eps=1e-11,iter_max=1000,toler_chol=coxph.control()$toler.chol,status='PASS')}))
stopifnot(all(fit_audit$minimum_covariance_eigenvalue>0))
write_tsv(fit_audit,'task5_fit_diagnostics.tsv','audit')
write_tsv(do.call(rbind,nested),'task5_nested_model_audit.tsv','audit')
write_tsv(effects,'10_cohort_primary_effects.tsv');write_tsv(meta_out,'11_random_effects_meta_analysis.tsv');write_tsv(interactions,'12_exploratory_interactions_bh.tsv');write_tsv(diagnostics,'13_proportional_hazards_diagnostics.tsv');write_tsv(do.call(rbind,feasibility),'14_interaction_feasibility.tsv');write_tsv(do.call(rbind,slopes),'15_subgroup_score_slopes.tsv')
write_tsv(do.call(rbind,checks),'task5_numeric_agreement.tsv','audit')
write_tsv(data.frame(tau2_closed_form=tau_independent,pooled_log_hr=mu,q=q,q_modified=qm,se=sem,lower95_log=ind_ci[1],upper95_log=ind_ci[2],prediction_lower95_log=ind_pi[1],prediction_upper95_log=ind_pi[2],status='PASS'),'task5_independent_meta_check.tsv','audit')
# Do not serialize coxph model/x/y objects: they retain patient-level data.
capture.output(sessionInfo(),file=file.path(out,'logs/task5_sessionInfo.txt'))
status<-if(any(interactions$status!='ESTIMATED'))'COMPLETE_WITH_DOCUMENTED_LIMITATION' else 'COMPLETE'
writeLines(c(paste0('status=',status),'original_gene_models_fitted=0','frozen_score_associations_only=TRUE','random_null_completion_dependency=FALSE','source_objects_serialized=FALSE'),file.path(out,'audit/task5_execution_status.txt'))
print(effects);print(meta_out);print(interactions);cat('\nTask5:',status,'\n')
