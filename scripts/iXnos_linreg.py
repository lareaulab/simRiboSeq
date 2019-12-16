import iXnos.interface as inter

### specify expt
base_dir = "/mnt/lareaulab/amok/iXnos/expts/"
expt = "scer_biasDelta_withBias_extraBase"
expt_dir = base_dir + expt
sam_fname = expt_dir + "/process/" + "scer_biasDelta_withBias_extraBase.transcript.mapped.wts.sam"
tr_codons_fname = expt_dir + "/process/" + "te_set_bounds.size.27.31.trunc.20.20.min_cts.200.min_cod.100.top.500.txt"
te_codons_fname = expt_dir + "/process/" + "te_set_bounds.size.27.31.trunc.20.20.min_cts.200.min_cod.100.top.500.txt"
outputs_fname = expt_dir + "/process/" + "outputs.size.27.31.txt"

### specify references
genome_data_dir = "/mnt/lareaulab/amok/iXnos/genome_data/"
gene_len_fname = genome_data_dir + "scer.transcripts.20cds20.lengths.txt"
gene_seq_fname = genome_data_dir + "scer.transcripts.20cds20.fa"

### specify model names
nb_model_names = [
        "cod_p0", "cod_p0_nt_p0p2", "cod_n3p2", "cod_n3p2_nt_n9p8",
        "cod_n5p4", "cod_n5p4_nt_n15p14", "cod_n7p5", "cod_n7p5_nt_n21p17"]
loo_model_names = ["nocod{0}_cod_n7p5_nt_n21p17".format(cod)
                    for cod in range(-7, 6)]

### compute neighborhood (nb) models:
for model_name in nb_model_names:
    print(model_name)
    model_params = model_name.split("_")
    cod_idxs = model_params[1]
    cod_start_sign = 1 if cod_idxs[0] == "p" else -1
    cod_start_idx = cod_start_sign * int(cod_idxs[1])
    cod_end_sign = 1 if cod_idxs[-2] == "p" else -1
    cod_end_idx = cod_end_sign * int(cod_idxs[-1])
    rel_cod_idxs = range(cod_start_idx, cod_end_idx + 1)
    if len(model_params) > 2:
        #Get start sign
        nt_idxs = model_params[3]
        nt_start_sign = 1 if nt_idxs[0] == "p" else -1
        nt_idxs = nt_idxs[1:]
        nt_idxs = nt_idxs.split("n")
        if len(nt_idxs) == 2:
            nt_end_sign = -1
        else:
            nt_idxs = nt_idxs[0]
            nt_idxs = nt_idxs.split("p")
            nt_end_sign = 1
        nt_start_idx = nt_start_sign * int(nt_idxs[0])
        nt_end_idx = nt_end_sign * int(nt_idxs[1])
        rel_nt_idxs = range(nt_start_idx, nt_end_idx + 1)
    else:
        rel_nt_idxs = []
    name = "lr_" + model_name
    wts, y_tr_hat, y_te_hat, y_tr, y_te = inter.make_linreg(
        expt_dir, name, gene_seq_fname, gene_len_fname, tr_codons_fname,
        te_codons_fname, outputs_fname, rel_cod_idxs=rel_cod_idxs,
        rel_nt_idxs=rel_nt_idxs)

### aggregate correlations
expt_linreg_dir = expt_dir + "/linreg"
result_dir = "/mnt/lareaulab/amok/iXnos/results/" + expt
out_fname = result_dir + "/feat_neighborhood_series/linreg_corrs.txt"
lr_nb_model_names = ["lr_" + model for model in nb_model_names]
inter.linreg.aggregate_linreg_corrs(expt_linreg_dir, lr_nb_model_names, out_fname)

### compute leave-one-out (loo) models
for model_name in loo_model_names:

    model_params = model_name.split("_")

    leaveout_idx = int(model_params[0][5:])

    rel_cod_idxs = range(-7, leaveout_idx) + range(leaveout_idx + 1, 6)
    rel_nt_idxs = range(-21, 3*leaveout_idx) + range(3*(leaveout_idx + 1), 18)

    name = "lr_" + model_name

    linreg = inter.make_linreg(
                expt_dir, name, gene_seq_fname, gene_len_fname, tr_codons_fname,
                te_codons_fname, outputs_fname, rel_cod_idxs=rel_cod_idxs, rel_nt_idxs=rel_nt_idxs,
                rel_struc_idxs=False, struc_fname=False)

### aggregate correlations
expt_linreg_dir = expt_dir + "/linreg"
result_dir = "/mnt/lareaulab/amok/iXnos/results/" + expt
out_fname = result_dir + "/leaveout_series/lr_leaveout_corrs.txt"
lr_loo_model_names = ["lr_" + model for model in loo_model_names]
inter.linreg.aggregate_linreg_corrs(expt_linreg_dir, lr_loo_model_names, out_fname)

