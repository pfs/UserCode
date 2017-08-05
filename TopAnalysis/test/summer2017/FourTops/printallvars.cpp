#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "TString.h"

void printallvars(const char *filename)
{
    if (!filename)
    {
        printf("No file specified. Exiting.\n");
    }
    TFile *file = new TFile(filename);
    TTree *tree = (TTree*) file->Get("TMVAanalysis");

    Float_t n_leptons, n_jets, n_bjets, n_non_bjets;
    Float_t jet_csv_1, jet_csv_2, jet_csv_3, jet_csv_4;
    Float_t n_mu_p, n_mu_m, n_ele_p, n_ele_m;
    Float_t angle_jets, angle_bjets;

    Float_t bdt, nn;

    tree->SetBranchAddress("n_leptons",&n_leptons);
    tree->SetBranchAddress("n_jets",&n_jets);
    tree->SetBranchAddress("n_bjets",&n_bjets);
    tree->SetBranchAddress("n_non_bjets",&n_non_bjets);

    tree->SetBranchAddress("jet_1_highest_csv",&jet_csv_1);
    tree->SetBranchAddress("jet_2_highest_csv",&jet_csv_2);
    tree->SetBranchAddress("jet_3_highest_csv",&jet_csv_3);
    tree->SetBranchAddress("jet_4_highest_csv",&jet_csv_4);

    tree->SetBranchAddress("n_mu_p",&n_mu_p);
    tree->SetBranchAddress("n_mu_m",&n_mu_m);
    tree->SetBranchAddress("n_ele_p",&n_ele_p);
    tree->SetBranchAddress("n_ele_m",&n_ele_m);

    tree->SetBranchAddress("jet_smallest_angle",&angle_jets);
    tree->SetBranchAddress("jet_smallest_angle_2b",&angle_bjets);

    tree->SetBranchAddress("bdt_discrim",&bdt);
    tree->SetBranchAddress("nn_discrim",&nn);

    FILE *out;
    out = fopen((string("out_")+string(filename)+string(".txt")).c_str(),"w");

    for (int i=0;i<tree->GetEntries();i++)
    {
        tree->GetEntry(i);
        fprintf(out, "%f\t", bdt);
        fprintf(out, "%f\t", nn);
        fprintf(out, "%f\t", n_leptons);
        fprintf(out, "%f\t", n_jets);
        fprintf(out, "%f\t", n_bjets);
        fprintf(out, "%f\t", n_non_bjets);
        fprintf(out, "%f\t", jet_csv_1);
        fprintf(out, "%f\t", jet_csv_2);
        fprintf(out, "%f\t", jet_csv_3);
        fprintf(out, "%f\t", jet_csv_4);
        fprintf(out, "%f\t", n_mu_p);
        fprintf(out, "%f\t", n_mu_m);
        fprintf(out, "%f\t", n_ele_p);
        fprintf(out, "%f\t", n_ele_m);
        fprintf(out, "%f\t", angle_jets);
        fprintf(out, "%f\t", angle_bjets);
        fprintf(out, "\n");
    }

    fclose(out);
}