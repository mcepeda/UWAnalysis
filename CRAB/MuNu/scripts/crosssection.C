double xsecTheo4F=138.9;  //46.3*3;
double xsecTheo5F=36703.2;  //12234.4*3

double DETECTORYIELD5F=905.4; // 1.16290925741195679e+03
double DETECTORYIELD4F=627.7; // 8.24458413273096085e+02
//double RFACTOR=1.330, ERROR_RFACTOR = 0.148; // Wbb4F 1m1e
//double RFACTOR=1.270, ERROR_RFACTOR = 0.154; // Wbb4F 3jet
double RFACTOR=1.351, ERROR_RFACTOR = 0.149; // Wbb5F 1m1e
//double RFACTOR=1.285, ERROR_RFACTOR = 0.149; // Wbb5F 3jet

double crosssection_4F=0.333233, error_4F=0.00149997;
double crosssection_5F=0.445021, error_5F=0.0148176;
double crosssection_5F_comb=0.472773, error_5F_comb=0.00700892;

double crosssection_4F_nofwd=0.364933, error_4F_nofwd=0.001569692;
double crosssection_5F_nofwd=0.492878, error_5F_nofwd=0.015594;
double crosssection_5F_comb_nofwd=0.524374, error_5F_comb_nofwd=0.0074137;

void crosssection(){

cout<<"==============="<<endl;
cout<<"A) WITH FWD VETO"<<endl;

cout<<" CROSSSECTION "<<endl;
cout<<"5F (with MPI): "<<crosssection_5F<<"+-"<<error_5F<<" pb"<<endl;
cout<<"5F_comb (with MPI): "<<crosssection_5F_comb<<"+-"<<error_5F_comb<<" pb"<<endl;
cout<<"4F (no MPI)  : "<<crosssection_4F<<"+-"<<error_4F<<" pb"<<endl;

cout<<"-----------"<<endl;

cout<<" ACCEPTANCE "<<endl;

double detector=DETECTORYIELD5F/19783.;   double error_detector=0;//40./19800;
double detector4F= DETECTORYIELD4F/19783.; double error_detector4F=0;//25./19800;

double acc5F=detector/crosssection_5F; 
double erroracc5F=acc5F*sqrt( (error_5F/crosssection_5F)*(error_5F/crosssection_5F) + (error_detector/detector)*(error_detector/detector));
double acc5F_comb=detector/crosssection_5F_comb; 
double erroracc5F_comb=acc5F_comb*sqrt( (error_5F_comb/crosssection_5F_comb)*(error_5F_comb/crosssection_5F_comb) + (error_detector/detector)*(error_detector/detector));
double acc4F=detector4F/crosssection_4F; 
double erroracc4F=acc4F*sqrt( (error_4F/crosssection_4F)*(error_4F/crosssection_4F) + (error_detector4F/detector4F)*(error_detector4F/detector4F));

cout<<"5F (with MPI): "<<acc5F*100<<"+-"<<erroracc5F*100<<" %"<<endl;
cout<<"5F_comb (with MPI): "<<acc5F_comb*100<<"+-"<<erroracc5F_comb*100<<" %"<<endl;
cout<<"4F (no MPI): "<<acc4F*100<<"+-"<<erroracc4F*100<<" %"<<endl;

double observed=DETECTORYIELD5F*RFACTOR; double error_observed=observed*sqrt( (40./DETECTORYIELD5F*40./DETECTORYIELD5F)+(ERROR_RFACTOR/RFACTOR)*(ERROR_RFACTOR/RFACTOR));
double observed_xsec=observed/acc5F_comb/19783;
double err_observed_xsec_exp=observed_xsec*error_observed/observed;
double err_observed_xsec_statacc=observed_xsec*erroracc5F_comb/acc5F_comb;
double err_observed_xsec_theoPDF=observed_xsec*0.01;
double err_observed_xsec_theoScale=observed_xsec*0.1;
double err_observed_xsec_theo=sqrt(err_observed_xsec_theoPDF*err_observed_xsec_theoPDF+err_observed_xsec_theoScale*err_observed_xsec_theoScale+err_observed_xsec_statacc*err_observed_xsec_statacc);

cout<<"-----------"<<endl;

cout<<"OBSERVED"<<endl;

cout<<observed_xsec<<"  (exp: "<<err_observed_xsec_exp<<"   theo: "<<err_observed_xsec_theo<<") pb"<<endl;

cout<<"==============="<<endl;
cout<<"B) WITHOUT FWD VETO"<<endl;

cout<<" CROSSSECTION "<<endl;
cout<<"5F_nofwd (with MPI): "<<crosssection_5F_nofwd<<"+-"<<error_5F_nofwd<<" pb"<<endl;
cout<<"5F_comb_nofwd (with MPI): "<<crosssection_5F_comb_nofwd<<"+-"<<error_5F_comb_nofwd<<" pb"<<endl;
cout<<"4F_nofwd (no MPI)  : "<<crosssection_4F_nofwd<<"+-"<<error_4F_nofwd<<" pb"<<endl;

cout<<"-----------"<<endl;

cout<<" ACCEPTANCE "<<endl;

double detector=DETECTORYIELD5F/19800;   double error_detector=0;//40./19800;
double detector4F= DETECTORYIELD4F/19800; double error_detector4F=0;//25./19800;

double acc5F_nofwd=detector/crosssection_5F_nofwd;
double erroracc5F_nofwd=acc5F_nofwd*sqrt( (error_5F_nofwd/crosssection_5F_nofwd)*(error_5F_nofwd/crosssection_5F_nofwd) + (error_detector/detector)*(error_detector/detector));
double acc5F_comb_nofwd=detector/crosssection_5F_comb_nofwd;
double erroracc5F_comb_nofwd=acc5F_comb_nofwd*sqrt( (error_5F_comb_nofwd/crosssection_5F_comb_nofwd)*(error_5F_comb_nofwd/crosssection_5F_comb_nofwd) + (error_detector/detector)*(error_detector/detector));
double acc4F_nofwd=detector4F/crosssection_4F_nofwd;
double erroracc4F_nofwd=acc4F_nofwd*sqrt( (error_4F_nofwd/crosssection_4F_nofwd)*(error_4F_nofwd/crosssection_4F_nofwd) + (error_detector4F/detector4F)*(error_detector4F/detector4F));

cout<<"5F (with MPI): "<<acc5F_nofwd*100<<"+-"<<erroracc5F_nofwd*100<<" %"<<endl;
cout<<"5F_comb (with MPI): "<<acc5F_comb_nofwd*100<<"+-"<<erroracc5F_comb_nofwd*100<<" %"<<endl;
cout<<"4F (no MPI): "<<acc4F_nofwd*100<<"+-"<<erroracc4F_nofwd*100<<" %"<<endl;

double observed=DETECTORYIELD5F*RFACTOR; double error_observed=observed*sqrt( (40./DETECTORYIELD5F*40./DETECTORYIELD5F)+(ERROR_RFACTOR/RFACTOR)*(ERROR_RFACTOR/RFACTOR));
double observed_xsec_nofwd=observed/acc5F_comb_nofwd/19800;
double err_observed_xsec_nofwd_exp=observed_xsec_nofwd*error_observed/observed;
double err_observed_xsec_nofwd_statacc=observed_xsec_nofwd*erroracc5F_comb/acc5F_comb;
double err_observed_xsec_nofwd_theoPDF=observed_xsec_nofwd*0.01;
double err_observed_xsec_nofwd_theoScale=observed_xsec_nofwd*0.1;
double err_observed_xsec_nofwd_theo=sqrt(err_observed_xsec_nofwd_theoPDF*err_observed_xsec_nofwd_theoPDF+err_observed_xsec_nofwd_theoScale*err_observed_xsec_nofwd_theoScale+err_observed_xsec_nofwd_statacc*err_observed_xsec_nofwd_statacc);

double observed_xsec4f_nofwd=observed/acc4F_nofwd/19800;
double err_observed_xsec4f_nofwd_exp=observed_xsec4f_nofwd*error_observed/observed;
double err_observed_xsec4f_nofwd_statacc=observed_xsec4f_nofwd*erroracc5F_comb/acc5F_comb;
double err_observed_xsec4f_nofwd_theoPDF=observed_xsec4f_nofwd*0.01;
double err_observed_xsec4f_nofwd_theoScale=observed_xsec4f_nofwd*0.1;
double err_observed_xsec4f_nofwd_theo=sqrt(err_observed_xsec4f_nofwd_theoPDF*err_observed_xsec4f_nofwd_theoPDF+err_observed_xsec4f_nofwd_theoScale*err_observed_xsec4f_nofwd_theoScale+err_observed_xsec4f_nofwd_statacc*err_observed_xsec4f_nofwd_statacc);



cout<<"-----------"<<endl;

cout<<"OBSERVED"<<endl;

cout<<"5F:  "<<observed_xsec_nofwd<<"  (exp: "<<err_observed_xsec_nofwd_exp<<"   theo: "<<err_observed_xsec_nofwd_theo<<") pb"<<endl;
cout<<"4F:  "<<observed_xsec4f_nofwd<<"  (exp: "<<err_observed_xsec4f_nofwd_exp<<"   theo: "<<err_observed_xsec4f_nofwd_theo<<") pb"<<endl;

cout<<"==============="<<endl;

}

