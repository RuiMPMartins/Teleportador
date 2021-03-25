#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH1.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TView.h"
#include "TGraphErrors.h"
#include <TGraph2DErrors.h>
#include <Spline3Interpolator.h>

#include "DataPoints.h"

using namespace std;


int main(int argc, char *argv[]){
    string line;
    double temp1;
    int Np = 24;
    ifstream infile("Files_Images_LET/Celula_Etanol/etanol01_18.dat", ios::in);

    //Coisas que se devem guardar
    vector<double> e1_t;
    vector<double> e1_I;
    vector<double> e1_V;

    for(int i=0; i<Np; i++){
    	getline(infile, line);
        stringstream ss;
        ss << line;
        string temp="";

        ss >> temp;
        if (stringstream(temp) >> temp1);
        e1_t.push_back(temp1);

        ss >> temp;
        if (stringstream(temp) >> temp1);
		e1_I.push_back(temp1);

        ss >> temp;
        if (stringstream(temp) >> temp1);
        e1_V.push_back(temp1);
    }

    infile.close();

    cout << "[FIM DA LEITURA DO FICHEIRO 1]" << endl;

    Np = 84;
    ifstream infile2("Files_Images_LET/Celula_Etanol/etanol02_18.dat", ios::in);

    //Coisas que se devem guardar
    vector<double> e2_t;
    vector<double> e2_I;
    vector<double> e2_V;

    for(int i=0; i<Np; i++){
        getline(infile2, line);
        stringstream ss;
        ss << line;
        string temp="";

        ss >> temp;
        if (stringstream(temp) >> temp1);
        e2_t.push_back(temp1);

        ss >> temp;
        if (stringstream(temp) >> temp1);
        e2_I.push_back(temp1);

        ss >> temp;
        if (stringstream(temp) >> temp1);
        e2_V.push_back(temp1);
    }

    infile2.close();

    cout << "[FIM DA LEITURA DO FICHEIRO 2]" << endl;

    Np = 48;
    ifstream infile3("Files_Images_LET/Celula_Etanol/etanol03_18.dat", ios::in);

    //Coisas que se devem guardar
    vector<double> e3_t;
    vector<double> e3_I;
    vector<double> e3_V;

    for(int i=0; i<Np; i++){
        getline(infile3, line);
        stringstream ss;
        ss << line;
        string temp="";

        ss >> temp;
        if (stringstream(temp) >> temp1);
        e3_t.push_back(temp1);

        ss >> temp;
        if (stringstream(temp) >> temp1);
        e3_I.push_back(temp1);

        ss >> temp;
        if (stringstream(temp) >> temp1);
        e3_V.push_back(temp1);
    }

    infile3.close();

    cout << "[FIM DA LEITURA DO FICHEIRO 3]" << endl;

    Np = 24;
    ifstream infile4("Files_Images_LET/Celula_Etanol/etanol04_18.dat", ios::in);

    //Coisas que se devem guardar
    vector<double> e4_t;
    vector<double> e4_I;
    vector<double> e4_V;

    for(int i=0; i<Np; i++){
        getline(infile4, line);
        stringstream ss;
        ss << line;
        string temp="";

        ss >> temp;
        if (stringstream(temp) >> temp1);
        e4_t.push_back(temp1);

        ss >> temp;
        if (stringstream(temp) >> temp1);
        e4_I.push_back(temp1);

        ss >> temp;
        if (stringstream(temp) >> temp1);
        e4_V.push_back(temp1);
    }

    infile4.close();

    cout << "[FIM DA LEITURA DO FICHEIRO 4]" << endl;

    Np = 43;
    double temp2;
    ifstream infile5("Files_Images_LET/Celula_Etanol/Fc_ph.txt", ios::in);

    //Coisas que se devem guardar
    vector<pair<double,double>> ph_fc;
    pair<double,double> temp_pair;

    for(int i=0; i<Np; i++){
        getline(infile5, line);
        stringstream ss;
        ss << line;
        string temp="";

        ss >> temp;
        if (stringstream(temp) >> temp1);

        ss >> temp;
        if (stringstream(temp) >> temp2);
        
        temp_pair = make_pair(temp2,temp1);
        ph_fc.push_back(temp_pair);
    }

    infile5.close();

    cout << "[FIM DA LEITURA DO FICHEIRO Fc-pH]" << endl;


    //Determinação do fator de conversão
    double e1_ph = 3.2, e2_ph = 3.31, e3_ph = 3.05, e4_ph = 3.08;
    double e1_fc, e2_fc, e3_fc, e4_fc;
    double e1_erfc, e2_erfc, e3_erfc, e4_erfc;

    Spline3Interpolator S3(ph_fc);

    e1_fc = S3.Interpolate(e1_ph);
    e2_fc = S3.Interpolate(e2_ph);
    e3_fc = S3.Interpolate(e3_ph);
    e4_fc = S3.Interpolate(e4_ph);

    S3.Draw();

    //erro Fc
    e1_erfc = 0.01*fabs(S3.Interpolate(e1_ph+1e-4)-S3.Interpolate(e1_ph-1e-4))/2e-4;
    e2_erfc = 0.01*fabs(S3.Interpolate(e2_ph+1e-4)-S3.Interpolate(e2_ph-1e-4))/2e-4;
    e3_erfc = 0.01*fabs(S3.Interpolate(e3_ph+1e-4)-S3.Interpolate(e3_ph-1e-4))/2e-4;
    e4_erfc = 0.01*fabs(S3.Interpolate(e4_ph+1e-4)-S3.Interpolate(e4_ph-1e-4))/2e-4;

    cout << endl;
    cout << "Fator de conversão Ensaio 1 = " << e1_fc << " +/- " << e1_erfc << " ( " << e1_erfc*100/e1_fc << " % )"<< endl;
    cout << "Fator de conversão Ensaio 2 = " << e2_fc << " +/- " << e2_erfc << " ( " << e2_erfc*100/e2_fc << " % )"<< endl;
    cout << "Fator de conversão Ensaio 3 = " << e3_fc << " +/- " << e3_erfc << " ( " << e3_erfc*100/e3_fc << " % )"<< endl;
    cout << "Fator de conversão Ensaio 4 = " << e4_fc << " +/- " << e4_erfc << " ( " << e4_erfc*100/e4_fc << " % )"<< endl;
    ////

    //Dados adicionais
    double e1_T = 292.65, e2_T = e1_T, e3_T3 = 313.15, e4_T = 328.15;      //Temperatura       (K)
    double e1_dt = 7923, e2_dt = 2479, e3_dt = 6416, e4_dt = 3394;             //Tempo total       (s)
    double e1_m = 12.5868, e2_m = 4.1842, e3_m = 11.7198, e4_m = 5.1886;   // Massa de efluente(kg)
    double e2_R = 60; //Resistencia constante no ensaior 2 (Ohm)

    //Menores divisoes da escala
    double er_m = 0.0001;

    //P_Eq e o seu erro
    double e1_P_Eq = 862.3*e1_fc*e1_m/((1+0.056*e1_fc)*e1_dt);
    double e2_P_Eq = 862.3*e2_fc*e2_m/((1+0.056*e2_fc)*e2_dt);
    double e3_P_Eq = 862.3*e3_fc*e3_m/((1+0.056*e3_fc)*e3_dt);
    double e4_P_Eq = 862.3*e4_fc*e4_m/((1+0.056*e4_fc)*e4_dt);

    double e1_er_P_Eq = fabs(26946875*e1_m*e1_erfc/(2*e1_dt*(7*e1_fc+125)*(7*e1_fc+125)))+fabs(215575*e1_fc*er_m/((14+250)*e1_dt));
    double e2_er_P_Eq = fabs(26946875*e2_m*e2_erfc/(2*e2_dt*(7*e2_fc+125)*(7*e2_fc+125)))+fabs(215575*e2_fc*er_m/((14+250)*e2_dt));
    double e3_er_P_Eq = fabs(26946875*e3_m*e3_erfc/(2*e3_dt*(7*e3_fc+125)*(7*e3_fc+125)))+fabs(215575*e3_fc*er_m/((14+250)*e3_dt));
    double e4_er_P_Eq = fabs(26946875*e4_m*e4_erfc/(2*e4_dt*(7*e4_fc+125)*(7*e4_fc+125)))+fabs(215575*e4_fc*er_m/((14+250)*e4_dt));

    cout << endl;
    cout << "P_Eq Ensaio 1 = " << e1_P_Eq << " +/- " << e1_er_P_Eq << " ( " << e1_er_P_Eq*100/e1_P_Eq << " % )"<< endl;
    cout << "P_Eq Ensaio 2 = " << e2_P_Eq << " +/- " << e2_er_P_Eq << " ( " << e2_er_P_Eq*100/e2_P_Eq << " % )"<< endl;
    cout << "P_Eq Ensaio 3 = " << e3_P_Eq << " +/- " << e3_er_P_Eq << " ( " << e3_er_P_Eq*100/e3_P_Eq << " % )"<< endl;
    cout << "P_Eq Ensaio 4 = " << e4_P_Eq << " +/- " << e4_er_P_Eq << " ( " << e4_er_P_Eq*100/e4_P_Eq << " % )"<< endl;
    cout << endl;
    ////


    TGraphErrors* gr_en1_V_I = new TGraphErrors();
    TGraphErrors* gr_en2_V_I = new TGraphErrors();
    TGraphErrors* gr_en3_V_I = new TGraphErrors();
    TGraphErrors* gr_en4_V_I = new TGraphErrors();

    TGraphErrors* gr_en1_P_t = new TGraphErrors();
    TGraphErrors* gr_en2_P_t = new TGraphErrors();
    TGraphErrors* gr_en3_P_t = new TGraphErrors();
    TGraphErrors* gr_en4_P_t = new TGraphErrors();

    TGraphErrors* gr_en1_R_t = new TGraphErrors();
    TGraphErrors* gr_en2_R_t = new TGraphErrors();
    TGraphErrors* gr_en3_R_t = new TGraphErrors();
    TGraphErrors* gr_en4_R_t = new TGraphErrors();

    TGraphErrors* gr_en1_V_t = new TGraphErrors();
    TGraphErrors* gr_en2_V_t = new TGraphErrors();
    TGraphErrors* gr_en3_V_t = new TGraphErrors();
    TGraphErrors* gr_en4_V_t = new TGraphErrors();

    TGraphErrors* gr_en1_I_t = new TGraphErrors();
    TGraphErrors* gr_en2_I_t = new TGraphErrors();
    TGraphErrors* gr_en3_I_t = new TGraphErrors();
    TGraphErrors* gr_en4_I_t = new TGraphErrors();

    TGraphErrors* gr_en1_P_R = new TGraphErrors();
    TGraphErrors* gr_en2_P_R = new TGraphErrors();
    TGraphErrors* gr_en3_P_R = new TGraphErrors();
    TGraphErrors* gr_en4_P_R = new TGraphErrors();

    TGraphErrors* gr_en1_P_R_S3 = new TGraphErrors();
    TGraphErrors* gr_en2_P_R_S3 = new TGraphErrors();
    TGraphErrors* gr_en3_P_R_S3 = new TGraphErrors();
    TGraphErrors* gr_en4_P_R_S3 = new TGraphErrors();

    TCanvas* c_en1_V_I = new TCanvas("c_en1_V_I", "c_en1_V_I", 800, 600);
    TCanvas* c_en2_V_I = new TCanvas("c_en2_V_I", "c_en2_V_I", 800, 600);
    TCanvas* c_en3_V_I = new TCanvas("c_en3_V_I", "c_en3_V_I", 800, 600);
    TCanvas* c_en4_V_I = new TCanvas("c_en4_V_I", "c_en4_V_I", 800, 600);
    TCanvas* c_V_I = new TCanvas("c_V_I", "c_V_I", 800, 600);

    TCanvas* c_en1_P_t = new TCanvas("c_en1_P_t", "c_en1_P_t", 800, 600);
    TCanvas* c_en2_P_t = new TCanvas("c_en2_P_t", "c_en2_P_t", 800, 600);
    TCanvas* c_en3_P_t = new TCanvas("c_en3_P_t", "c_en3_P_t", 800, 600);
    TCanvas* c_en4_P_t = new TCanvas("c_en4_P_t", "c_en4_P_t", 800, 600);

    TCanvas* c_en1_R_t = new TCanvas("c_en1_R_t", "c_en1_R_t", 800, 600);
    TCanvas* c_en2_R_t = new TCanvas("c_en2_R_t", "c_en2_R_t", 800, 600);
    TCanvas* c_en3_R_t = new TCanvas("c_en3_R_t", "c_en3_R_t", 800, 600);
    TCanvas* c_en4_R_t = new TCanvas("c_en4_R_t", "c_en4_R_t", 800, 600);

    TCanvas* c_en1_V_t = new TCanvas("c_en1_V_t", "c_en1_V_t", 800, 600);
    TCanvas* c_en2_V_t = new TCanvas("c_en2_V_t", "c_en2_V_t", 800, 600);
    TCanvas* c_en3_V_t = new TCanvas("c_en3_V_t", "c_en3_V_t", 800, 600);
    TCanvas* c_en4_V_t = new TCanvas("c_en4_V_t", "c_en4_V_t", 800, 600);

    TCanvas* c_en1_I_t = new TCanvas("c_en1_I_t", "c_en1_I_t", 800, 600);
    TCanvas* c_en2_I_t = new TCanvas("c_en2_I_t", "c_en2_I_t", 800, 600);
    TCanvas* c_en3_I_t = new TCanvas("c_en3_I_t", "c_en3_I_t", 800, 600);
    TCanvas* c_en4_I_t = new TCanvas("c_en4_I_t", "c_en4_I_t", 800, 600);

    TCanvas* c_en1_P_R = new TCanvas("c_en1_P_R", "c_en1_P_R", 800, 600);
    TCanvas* c_en2_P_R = new TCanvas("c_en2_P_R", "c_en2_P_R", 800, 600);
    TCanvas* c_en3_P_R = new TCanvas("c_en3_P_R", "c_en3_P_R", 800, 600);
    TCanvas* c_en4_P_R = new TCanvas("c_en4_P_R", "c_en4_P_R", 800, 600);

    vector<pair<double,double>> e1_P_R;
    vector<pair<double,double>> e3_P_R;
    vector<pair<double,double>> e4_P_R;

    for(int i=0; i<24; i++){
        gr_en1_V_I->SetPoint(i,e1_I[i],e1_V[i]);
        gr_en1_V_I->SetPointError(i,e1_I[i]*2e-3,e1_V[i]*2e-3);

        gr_en1_P_t->SetPoint(i,e1_t[i],e1_I[i]*e1_V[i]);
        gr_en1_P_t->SetPointError(i,0,e1_V[i]*e1_I[i]*2e-3+e1_I[i]*e1_V[i]*2e-3);

        gr_en1_R_t->SetPoint(i,e1_t[i],e1_V[i]/e1_I[i]);
        gr_en1_R_t->SetPointError(i,0,e1_V[i]*2e-3/e1_I[i]+e1_V[i]*e1_I[i]*2e-3/(e1_I[i]*e1_I[i]));

        gr_en1_V_t->SetPoint(i,e1_t[i],e1_V[i]);
        gr_en1_V_t->SetPointError(i,0,e1_V[i]*2e-3);

        gr_en1_I_t->SetPoint(i,e1_t[i],e1_I[i]);
        gr_en1_I_t->SetPointError(i,0,e1_I[i]*2e-3);

        gr_en1_P_R->SetPoint(i,e1_V[i]/e1_I[i],e1_I[i]*e1_V[i]);
        gr_en1_P_R->SetPointError(i,e1_V[i]*2e-3/e1_I[i]+e1_V[i]*e1_I[i]*2e-3/(e1_I[i]*e1_I[i]),e1_V[i]*e1_I[i]*2e-3+e1_I[i]*e1_V[i]*2e-3);

        e1_P_R.push_back(make_pair(e1_V[i]/e1_I[i],e1_I[i]*e1_V[i]));
    }
    Spline3Interpolator S3_e1_P_R(e1_P_R);

    for(int i=0; i<84; i++){
        gr_en2_V_I->SetPoint(i,e2_I[i],e2_V[i]);
        gr_en2_V_I->SetPointError(i,e2_I[i]*2e-3,e2_V[i]*2e-3);

        gr_en2_P_t->SetPoint(i,e2_t[i],e2_I[i]*e2_V[i]);
        gr_en2_P_t->SetPointError(i,0,e2_V[i]*e2_I[i]*2e-3+e2_I[i]*e2_V[i]*2e-3);
        
        gr_en2_R_t->SetPoint(i,e2_t[i],e2_V[i]/e2_I[i]);
        gr_en2_R_t->SetPointError(i,0,e2_V[i]*2e-3/e2_I[i]+e2_V[i]*e2_I[i]*2e-3/(e2_I[i]*e2_I[i]));

        gr_en2_V_t->SetPoint(i,e2_t[i],e2_V[i]);
        gr_en2_V_t->SetPointError(i,0,e2_V[i]*2e-3);

        gr_en2_I_t->SetPoint(i,e2_t[i],e2_I[i]);
        gr_en2_I_t->SetPointError(i,0,e2_I[i]*2e-3);

        gr_en2_P_R->SetPoint(i,e2_V[i]/e2_I[i],e2_I[i]*e2_V[i]);
        gr_en2_P_R->SetPointError(i,e2_V[i]*2e-3/e2_I[i]+e2_V[i]*e2_I[i]*2e-3/(e2_I[i]*e2_I[i]),e2_V[i]*e2_I[i]*2e-3+e2_I[i]*e2_V[i]*2e-3);
    }

    for(int i=0; i<48; i++){
        gr_en3_V_I->SetPoint(i,e3_I[i],e3_V[i]);
        gr_en3_V_I->SetPointError(i,e3_I[i]*2e-3,e3_V[i]*2e-3);

        gr_en3_P_t->SetPoint(i,e3_t[i],e3_I[i]*e3_V[i]);
        gr_en3_P_t->SetPointError(i,0,e3_V[i]*e3_I[i]*2e-3+e3_I[i]*e3_V[i]*2e-3);
    
        gr_en3_R_t->SetPoint(i,e3_t[i],e3_V[i]/e3_I[i]);
        gr_en3_R_t->SetPointError(i,0,e3_V[i]*2e-3/e3_I[i]+e3_V[i]*e3_I[i]*2e-3/(e3_I[i]*e3_I[i]));

        gr_en3_V_t->SetPoint(i,e3_t[i],e3_V[i]);
        gr_en3_V_t->SetPointError(i,0,e3_V[i]*2e-3);

        gr_en3_I_t->SetPoint(i,e3_t[i],e3_I[i]);
        gr_en3_I_t->SetPointError(i,0,e3_I[i]*2e-3);

        gr_en3_P_R->SetPoint(i,e3_V[i]/e3_I[i],e3_I[i]*e3_V[i]);
        gr_en3_P_R->SetPointError(i,e3_V[i]*2e-3/e3_I[i]+e3_V[i]*e3_I[i]*2e-3/(e3_I[i]*e3_I[i]),e3_V[i]*e3_I[i]*2e-3+e3_I[i]*e3_V[i]*2e-3);
        
        e3_P_R.push_back(make_pair(e3_V[i]/e3_I[i],e3_I[i]*e3_V[i]));
    }
    Spline3Interpolator S3_e3_P_R(e3_P_R);

    for(int i=0; i<24; i++){
        gr_en4_V_I->SetPoint(i,e4_I[i],e4_V[i]);
        gr_en4_V_I->SetPointError(i,e4_I[i]*2e-3,e4_V[i]*2e-3);

        gr_en4_P_t->SetPoint(i,e4_t[i],e4_I[i]*e4_V[i]);
        gr_en4_P_t->SetPointError(i,0,e4_V[i]*e4_I[i]*2e-3+e4_I[i]*e4_V[i]*2e-3);
        
        gr_en4_R_t->SetPoint(i,e4_t[i],e4_V[i]/e4_I[i]);
        gr_en4_R_t->SetPointError(i,0,e4_V[i]*2e-3/e4_I[i]+e4_V[i]*e4_I[i]*2e-3/(e4_I[i]*e4_I[i]));

        gr_en4_V_t->SetPoint(i,e4_t[i],e4_V[i]);
        gr_en4_V_t->SetPointError(i,0,e4_V[i]*2e-3);

        gr_en4_I_t->SetPoint(i,e4_t[i],e4_I[i]);
        gr_en4_I_t->SetPointError(i,0,e4_I[i]*2e-3);

        gr_en4_P_R->SetPoint(i,e4_V[i]/e4_I[i],e4_I[i]*e4_V[i]);
        gr_en4_P_R->SetPointError(i,e4_V[i]*2e-3/e4_I[i]+e4_V[i]*e4_I[i]*2e-3/(e4_I[i]*e4_I[i]),e4_V[i]*e4_I[i]*2e-3+e4_I[i]*e4_V[i]*2e-3);
    
        e4_P_R.push_back(make_pair(e4_V[i]/e4_I[i],e4_I[i]*e4_V[i]));
    }
    Spline3Interpolator S3_e4_P_R(e4_P_R);

    TF1 *f = new TF1("f",[&](double*x, double *p){ 
        return p[0]-p[1]*x[0]; }, 0, 300,2);
    f->SetLineColor(kRed-3);

    TF1 *f_1 = new TF1("f_1",[&](double*x, double *p){ 
        return p[0]-p[1]*x[0]; }, 0, 300,2);

    TF1 *f_3 = new TF1("f_3",[&](double*x, double *p){ 
        return p[0]-p[1]*x[0]; }, 0, 300,2);

    TF1 *f_4 = new TF1("f_4",[&](double*x, double *p){ 
        return p[0]-p[1]*x[0]; }, 0, 300,2);

    c_en1_V_I->cd();
    gr_en1_V_I->SetLineColor(13);
    gr_en1_V_I->Fit(f_1,"fit 0");
    gr_en1_V_I->GetYaxis()->SetTitleOffset(1.4);
    gr_en1_V_I->SetTitle(";Corrente El#acute{e}trica (A); Tens#tilde{a}o El#acute{e}trica (V)");
    gr_en1_V_I->Draw("AP ERR");
    c_en1_V_I->SaveAs("Files_Images_LET/Celula_Etanol/En1_V_I.pdf");

    c_en2_V_I->cd();
    gr_en2_V_I->SetLineColor(13);
    gr_en2_V_I->Fit(f,"fit 0");
    TGraph* gr_fit_V_I_2 = new TGraph;
    for(int i=0; i<120; i++)
        gr_fit_V_I_2->SetPoint(i,0.0030+i*(0.0006/100),f->Eval(0.0030+i*(0.0006/100)));
    gr_en2_V_I->GetYaxis()->SetTitleOffset(1.5);
    gr_en2_V_I->SetTitle(";Corrente El#acute{e}trica (A); Tens#tilde{a}o El#acute{e}trica (V)");
    gr_en2_V_I->SetLineColor(kTeal+3);
    gr_fit_V_I_2->SetLineColor(kAzure-1);
    gr_en2_V_I->Draw("AP ERR");
    gr_fit_V_I_2->Draw("L same");
    c_en2_V_I->SaveAs("Files_Images_LET/Celula_Etanol/En2_V_I.pdf");

    c_en3_V_I->cd();
    gr_en3_V_I->SetLineColor(13);
    gr_en3_V_I->Fit(f_3,"fit 0");
    gr_en3_V_I->GetYaxis()->SetTitleOffset(1.4);
    gr_en3_V_I->SetTitle(";Corrente El#acute{e}trica (A); Tens#tilde{a}o El#acute{e}trica (V)");
    gr_en3_V_I->Draw("AP ERR");
    c_en3_V_I->SaveAs("Files_Images_LET/Celula_Etanol/En3_V_I.pdf");

    c_en4_V_I->cd();
    gr_en4_V_I->SetLineColor(13);
    gr_en4_V_I->Fit(f_4,"fit 0");
    gr_en4_V_I->GetYaxis()->SetTitleOffset(1.4);
    gr_en4_V_I->SetTitle(";Corrente El#acute{e}trica (A); Tens#tilde{a}o El#acute{e}trica (V)");
    gr_en4_V_I->Draw("AP ERR");
    gr_en3_V_I->Draw("P ERR same");
    gr_en1_V_I->Draw("P ERR same");
    c_en4_V_I->SaveAs("Files_Images_LET/Celula_Etanol/En4_V_I.pdf");

    c_V_I->cd();

    TGraph* gr_fit_V_I_1 = new TGraph;
    TGraph* gr_fit_V_I_3 = new TGraph;
    TGraph* gr_fit_V_I_4 = new TGraph;

    for(int i=0; i<120; i++){
        gr_fit_V_I_1->SetPoint(i,i*(0.008/100),f_1->Eval(i*(0.008/100)));
        gr_fit_V_I_3->SetPoint(i,i*(0.008/100),f_3->Eval(i*(0.008/100)));
        gr_fit_V_I_4->SetPoint(i,i*(0.008/100),f_4->Eval(i*(0.008/100)));
    }
    gr_fit_V_I_1->SetLineColor(kTeal-3);
    gr_fit_V_I_3->SetLineColor(kViolet-3);
    gr_fit_V_I_4->SetLineColor(kRed-3);
    gr_en1_V_I->SetLineColor(kTeal+3);
    gr_en3_V_I->SetLineColor(kViolet+3);
    gr_en4_V_I->SetLineColor(kRed+3);

    TGraph* grQuadra = new TGraph();
    grQuadra->SetPoint(0,0.008,0.15);
    grQuadra->SetPoint(1,0.008,0.5);
    grQuadra->SetPoint(2,0,0.15);
    grQuadra->SetPoint(3,0,0.5);
    grQuadra->SetMarkerColor(kWhite);
    grQuadra->SetTitle(";Corrente El#acute{e}trica (A); Tens#tilde{a}o El#acute{e}trica (V)");

    grQuadra->Draw("AP");
    gr_en4_V_I->Draw("P ERR same");
    gr_en3_V_I->Draw("P ERR same");
    gr_en1_V_I->Draw("P ERR same");
    gr_fit_V_I_1->Draw("L same");
    gr_fit_V_I_3->Draw("L same");
    gr_fit_V_I_4->Draw("L same");

    c_V_I->SaveAs("Files_Images_LET/Celula_Etanol/V_I.pdf");
    ////


    c_en1_P_t->cd();
    gr_en1_P_t->GetYaxis()->SetTitleOffset(1.55);
    gr_en1_P_t->SetTitle(";Tempo (s); Pot#hat{e}ncia (W)");
    gr_en1_P_t->Draw("AP ERR");
    c_en1_P_t->SaveAs("Files_Images_LET/Celula_Etanol/En1_P_t.pdf");

    c_en2_P_t->cd();
    gr_en2_P_t->GetYaxis()->SetTitleOffset(1.55);
    gr_en2_P_t->SetTitle(";Tempo (s); Pot#hat{e}ncia (W)");
    gr_en2_P_t->Draw("AP ERR");
    c_en2_P_t->SaveAs("Files_Images_LET/Celula_Etanol/En2_P_t.pdf");

    c_en3_P_t->cd();
    gr_en3_P_t->GetYaxis()->SetTitleOffset(1.55);
    gr_en3_P_t->SetTitle(";Tempo (s); Pot#hat{e}ncia (W)");
    gr_en3_P_t->Draw("AP ERR");
    c_en3_P_t->SaveAs("Files_Images_LET/Celula_Etanol/En3_P_t.pdf");

    c_en4_P_t->cd();
    gr_en4_P_t->GetYaxis()->SetTitleOffset(1.55);
    gr_en4_P_t->SetTitle(";Tempo (s); Pot#hat{e}ncia (W)");
    gr_en4_P_t->Draw("AP ERR");
    c_en4_P_t->SaveAs("Files_Images_LET/Celula_Etanol/En4_P_t.pdf");


    c_en1_R_t->cd();
    gr_en1_R_t->GetYaxis()->SetTitleOffset(1.4);
    gr_en1_R_t->SetTitle(";Tempo (s); Resist#hat{e}ncia (#Omega)");
    gr_en1_R_t->Draw("AP ERR");
    c_en1_R_t->SaveAs("Files_Images_LET/Celula_Etanol/En1_R_t.pdf");

    c_en2_R_t->cd();
    gr_en2_R_t->GetYaxis()->SetTitleOffset(1.4);
    gr_en2_R_t->SetTitle(";Tempo (s); Resist#hat{e}ncia (#Omega)");
    gr_en2_R_t->Draw("AP ERR");
    c_en2_R_t->SaveAs("Files_Images_LET/Celula_Etanol/En2_R_t.pdf");

    c_en3_R_t->cd();
    gr_en3_R_t->GetYaxis()->SetTitleOffset(1.4);
    gr_en3_R_t->SetTitle(";Tempo (s); Resist#hat{e}ncia (#Omega)");
    gr_en3_R_t->Draw("AP ERR");
    c_en3_R_t->SaveAs("Files_Images_LET/Celula_Etanol/En3_R_t.pdf");

    c_en4_R_t->cd();
    gr_en4_R_t->GetYaxis()->SetTitleOffset(1.4);
    gr_en4_R_t->SetTitle(";Tempo (s); Resist#hat{e}ncia (#Omega)");
    gr_en4_R_t->Draw("AP ERR");
    c_en4_R_t->SaveAs("Files_Images_LET/Celula_Etanol/En4_R_t.pdf");


    TF1 *f3 = new TF1("f3",[&](double*x, double *p){ 
        return p[0]+ p[1]*exp(-p[2]*x[0]); }, 0, 300,3);
    f3->SetLineColor(kRed-3);

    double R_media = 0;
    TGraphErrors* gr_en2_V_t_fit = new TGraphErrors();
    for(int i=0; i<84; i++){
        if(i>=4){
            gr_en2_V_t_fit->SetPoint(i-4,e2_t[i],e2_V[i]);
            gr_en2_V_t_fit->SetPointError(i-4,0,e2_V[i]*2e-3);
        }

        R_media += e2_V[i]/e2_I[i];
    }
    R_media /= 84;

    double e2_er_R = 0;
    for(int i=0; i<84; i++)
        e2_er_R += (R_media - e2_V[i]/e2_I[i])*(R_media - e2_V[i]/e2_I[i]);
    e2_er_R = sqrt(e2_er_R/83);

    gr_en2_V_t_fit->Fit(f3,"fit");

    double en2_V_eq = f3->GetParameter(0);
    double alpha = f3->GetParameter(1);
    double beta = f3->GetParameter(2);

    double I_eq = en2_V_eq/R_media;
    double er_I_eq = f3->GetParError(0)/R_media+en2_V_eq*e2_er_R/(R_media*R_media);
    double P_eq = en2_V_eq*I_eq;
    double er_P_eq = I_eq*f3->GetParError(0)+en2_V_eq*er_I_eq;

    cout << endl;
    cout << "Para a Temperatura de 19.5 °C, na Resistencia de " << R_media << " +/- " << e2_er_R << " ( " << e2_er_R*100/R_media << " % )"<< endl;
    cout << "A Voltagem de equilibrio é " << en2_V_eq << " +/- " << f3->GetParError(0) << " ( " << f3->GetParError(0)*100/en2_V_eq << " % )"<< endl;
    cout << "Isto implica que neste ensaio I_eq é " << I_eq << " +/- " << er_I_eq << " ( " << er_I_eq*100/I_eq << " % )"<< endl;
    cout << "A Potência de equilibrio é " << P_eq << " +/- " << er_P_eq << " ( " << er_P_eq*100/P_eq << " % )"<< endl;
    cout << "Assim o Rendimento de equilibrio é " << P_eq/e2_P_Eq << " +/- " << max(fabs(P_eq/e2_P_Eq-(P_eq-er_P_eq)/(e2_P_Eq+e2_er_P_Eq)),fabs(P_eq/e2_P_Eq-(P_eq+er_P_eq)/(e2_P_Eq-e2_er_P_Eq))) << " ( " << max(fabs(P_eq/e2_P_Eq-(P_eq-er_P_eq)/(e2_P_Eq+e2_er_P_Eq)),fabs(P_eq/e2_P_Eq-(P_eq+er_P_eq)/(e2_P_Eq-e2_er_P_Eq)))*100/(P_eq/e2_P_Eq) << " % )"<< endl;
    cout << endl;

    c_en1_V_t->cd();
    gr_en1_V_t->GetYaxis()->SetTitleOffset(1.4);
    gr_en1_V_t->SetTitle(";Tempo (s); Tens#tilde{a}o El#acute{e}trica (V)");
    gr_en1_V_t->Draw("AP ERR");
    gr_en1_I_t->Draw("P ERR same");
    c_en1_V_t->SaveAs("Files_Images_LET/Celula_Etanol/En1_V_t.pdf");

    c_en2_V_t->cd();
    gr_en2_V_t->GetYaxis()->SetTitleOffset(1.5);
    gr_en2_V_t->SetTitle(";Tempo (s); Tens#tilde{a}o El#acute{e}trica (V)");
    gr_en2_V_t->Draw("AP ERR");
    gr_en2_V_t_fit->Draw("P same");
    c_en2_V_t->SaveAs("Files_Images_LET/Celula_Etanol/En2_V_t.pdf");

    c_en3_V_t->cd();
    gr_en3_V_t->GetYaxis()->SetTitleOffset(1.4);
    gr_en3_V_t->SetTitle(";Tempo (s); Tens#tilde{a}o El#acute{e}trica (V)");
    gr_en3_V_t->Draw("AP ERR");
    c_en3_V_t->SaveAs("Files_Images_LET/Celula_Etanol/En3_V_t.pdf");

    c_en4_V_t->cd();
    gr_en4_V_t->GetYaxis()->SetTitleOffset(1.4);
    gr_en4_V_t->SetTitle(";Tempo (s); Tens#tilde{a}o El#acute{e}trica (V)");
    gr_en4_V_t->Draw("AP ERR");
    c_en4_V_t->SaveAs("Files_Images_LET/Celula_Etanol/En4_V_t.pdf");


    c_en1_I_t->cd();
    gr_en1_I_t->GetYaxis()->SetTitleOffset(1.4);
    gr_en1_I_t->SetTitle(";Tempo (s); Corrente El#acute{e}trica (A)");
    gr_en1_I_t->Draw("AP ERR");
    c_en1_I_t->SaveAs("Files_Images_LET/Celula_Etanol/En1_I_t.pdf");

    c_en2_I_t->cd();
    gr_en2_I_t->GetYaxis()->SetTitleOffset(1.4);
    gr_en2_I_t->SetTitle(";Tempo (s); Corrente El#acute{e}trica (A)");
    gr_en2_I_t->Draw("AP ERR");
    c_en2_I_t->SaveAs("Files_Images_LET/Celula_Etanol/En2_I_t.pdf");

    c_en3_I_t->cd();
    gr_en3_I_t->GetYaxis()->SetTitleOffset(1.4);
    gr_en3_I_t->SetTitle(";Tempo (s); Corrente El#acute{e}trica (A)");
    gr_en3_I_t->Draw("AP ERR");
    c_en3_I_t->SaveAs("Files_Images_LET/Celula_Etanol/En3_I_t.pdf");

    c_en4_I_t->cd();
    gr_en4_I_t->GetYaxis()->SetTitleOffset(1.4);
    gr_en4_I_t->SetTitle(";Tempo (s); Corrente El#acute{e}trica (A)");
    gr_en4_I_t->Draw("AP ERR");
    c_en4_I_t->SaveAs("Files_Images_LET/Celula_Etanol/En4_I_t.pdf");



    TF1 *f2 = new TF1("f",[&](double*x, double *p){ 
        return x[0]*p[0]*p[0]/((x[0]+p[1])*(x[0]+p[1])); }, 0, 300,2);
    f2->SetLineColor(kRed-3);
    f2->SetParameter(0,f->GetParameter(0));
    f2->SetParameter(1,f->GetParameter(1));

    for(int i=-20; i<220; i++){
        gr_en1_P_R_S3->SetPoint(i,e1_V[0]/e1_I[0]+i*(e1_V[23]/e1_I[23]-e1_V[0]/e1_I[0])/200,S3_e1_P_R.Interpolate(e1_V[0]/e1_I[0]+i*(e1_V[23]/e1_I[23]-e1_V[0]/e1_I[0])/200));
        gr_en3_P_R_S3->SetPoint(i,e3_V[0]/e3_I[0]+i*(e3_V[47]/e3_I[47]-e3_V[0]/e3_I[0])/200,S3_e3_P_R.Interpolate(e3_V[0]/e3_I[0]+i*(e3_V[47]/e3_I[47]-e3_V[0]/e3_I[0])/200));
        gr_en4_P_R_S3->SetPoint(i,e4_V[0]/e4_I[0]+i*(e4_V[23]/e4_I[23]-e4_V[0]/e4_I[0])/200,S3_e4_P_R.Interpolate(e4_V[0]/e4_I[0]+i*(e4_V[23]/e4_I[23]-e4_V[0]/e4_I[0])/200));
    }

    TF1 *f2_1 = new TF1("f",[&](double*x, double *p){ 
        return S3_e1_P_R.Interpolate(x[0]); }, 10, 300,0);
    TF1 *f2_3 = new TF1("f",[&](double*x, double *p){ 
        return S3_e3_P_R.Interpolate(x[0]); }, 10, 300,0);
    TF1 *f2_4 = new TF1("f",[&](double*x, double *p){ 
        return S3_e4_P_R.Interpolate(x[0]); }, 25, 300,0);

    c_en1_P_R->cd();
    gr_en1_P_R->SetLineColor(13);
    gr_en1_P_R->Fit(f2,"fit");
    gr_en1_P_R->GetYaxis()->SetTitleOffset(1.55);
    gr_en1_P_R->SetTitle(";Resist#hat{e}ncia (#Omega);  Pot#hat{e}ncia (W)");
    gr_en1_P_R->Draw("AP ERR");
    gr_en1_P_R_S3->SetLineColor(kGreen+3);
    gr_en1_P_R_S3->Draw("L same");
    c_en1_P_R->SaveAs("Files_Images_LET/Celula_Etanol/En1_P_R.pdf");

    cout << endl;
    cout << "A 19.5 °C a Resistencia ótima é " << f2->GetParameter(0) << " ( " << f2->GetMaximum() << " Watts) "  <<  " pelo fit e " << f2_1->GetX(f2_1->GetMaximum()-f2_1->GetMaximum()*0.0005) << " +/- 0.281956 ( " << f2_1->GetMaximum() << " +/- 4.34979e-06 Watts) "  << " pela interpolação Spline3" << endl;
    cout << endl;

    /*
    c_en2_P_R->cd();
    gr_en2_P_R->SetLineColor(13);
    gr_en2_P_R->Fit(f2,"fit");
    gr_en2_P_R->GetYaxis()->SetTitleOffset(1.55);
    gr_en2_P_R->SetTitle(";Resist#hat{e}ncia (#Omega);  Pot#hat{e}ncia (W)");
    gr_en2_P_R->Draw("AP ERR");
    c_en2_P_R->SaveAs("Files_Images_LET/Celula_Etanol/En2_P_R.pdf");*/

    c_en3_P_R->cd();
    gr_en3_P_R->SetLineColor(13);
    gr_en3_P_R->Fit(f2,"fit");
    gr_en3_P_R->GetYaxis()->SetTitleOffset(1.55);
    gr_en3_P_R->SetTitle(";Resist#hat{e}ncia (#Omega);  Pot#hat{e}ncia (W)");
    gr_en3_P_R->Draw("AP ERR");
    gr_en3_P_R_S3->SetLineColor(kGreen+3);
    gr_en3_P_R_S3->Draw("L same");
    c_en3_P_R->SaveAs("Files_Images_LET/Celula_Etanol/En3_P_R.pdf");

    cout << endl;
    cout << "A 40 °C a Resistencia ótima é " << f2->GetX(f2->GetMaximum()-f2->GetMaximum()*0.001) << " ( " << f2->GetMaximum() << " Watts) "  <<  " pelo fit e " << f2_3->GetX(f2_3->GetMaximum()-f2_3->GetMaximum()*0.001) << " +/- 0.260557 ( " << f2_3->GetMaximum() << " +/= 4.32327e-06 Watts) "  << " pela interpolação Spline3" << endl;
    cout << endl;

    c_en4_P_R->cd();
    gr_en4_P_R->SetLineColor(13);
    gr_en4_P_R->Fit(f2,"fit");
    gr_en4_P_R->GetYaxis()->SetTitleOffset(1.55);
    gr_en4_P_R->SetTitle(";Resist#hat{e}ncia (#Omega);  Pot#hat{e}ncia (W)");
    gr_en4_P_R->Draw("AP ERR");
    gr_en4_P_R_S3->SetLineColor(kGreen+3);
    gr_en4_P_R_S3->Draw("L same");
    c_en4_P_R->SaveAs("Files_Images_LET/Celula_Etanol/En4_P_R.pdf");
    ////

    cout << endl;
    cout << "A 55 °C a Resistencia ótima é " << f2->GetX(f2->GetMaximum()-f2->GetMaximum()*0.002) << " ( " << f2->GetMaximum() << " Watts) "  <<  " pelo fit e " << f2_4->GetX(f2_4->GetMaximum()-f2_4->GetMaximum()*0.0005) << " +/- 0.162099 ( " << f2_4->GetMaximum() << " +/- 1.18183e-05 Watts) "  << " pela interpolação Spline3" << endl;
    cout << endl;

    return 0;
}