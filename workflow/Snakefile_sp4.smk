include: "rules.smk"


metadata = pd.read_csv("metadata/runs_sp4.csv")


rule all_sp4:
    input:
        plots_per_run=expand( [
            "assets/sp4/plots/overview/{dataset}.pdf",
            "assets/sp4/plots/overview_zoom/{dataset}.pdf",
            "assets/sp4/plots/an_trajectories/{dataset}.pdf",
            "assets/sp4/plots/free_energy/free_energy_{dataset}.pdf",
            ], dataset=list_datasets(metadata) ),
        plots_Nt=expand( [
            "assets/sp4/plots/an_volume_Nt{Nt}.pdf",
            "assets/sp4/plots/entropy_Nt{Nt}.pdf",
            "assets/sp4/plots/critical_beta_ratios_Nt{Nt}.pdf",
            "assets/sp4/plots/plaquette_distribution_Nt{Nt}_volumes.pdf",
            "assets/sp4/plots/polyakov_susceptibility_Nt{Nt}.pdf", 
            "assets/sp4/plots/binder_cumulant_Nt{Nt}.pdf", 
            "assets/sp4/plots/specific_heat_Nt{Nt}.pdf", 
            ], Nt=[4, 5] ),
        double_gaussian=expand(
            "assets/sp4/plots/plaquette_distribution/pd_{dataset}_{ratio}.pdf",
            dataset=list_datasets(metadata),
            ratio=["1:1", "2:1", "1:2"] ),
        polyakov_loop_hist=expand( [
            "assets/sp4/plots/polyakov_loop_hist/polyakov_{dataset}.pdf",
            "assets/sp4/plots/polyakov_loop_compact/polyakov_{dataset}.pdf",
            ], dataset=list_datasets_poly(metadata), ),
        an_replicas=expand(
            "assets/sp4/plots/an_replicas_Nt{Nt}_Ns{Ns}.pdf", zip,
            Nt=[5, 5],
            Ns=[48, 56],
        ),
        tables = (
            "assets/sp4/tables/runs.tex",
            "assets/sp4/tables/beta_critical.tex",
        ),
        csv_files = (
            "data_assets/sp4/critical_beta_cumulants.csv",
            "data_assets/sp4/critical_beta_1:1.csv",
            "data_assets/sp4/critical_beta_2:1.csv",
            "data_assets/sp4/critical_beta_1:2.csv",
            "data_assets/sp4/supercooling_param.csv",
            "data_assets/sp4/scale_setting.csv",
            "data_assets/sp4/sp4_w0.csv",
        ),
        polyakov_loop_vs_beta="data_assets/sp4/polyakov/polyakov_loop_data.hdf5",
        surface_tension="assets/sp4/plots/surface_tension_term.pdf",
        wf_plot="assets/sp4/plots/vis_W_E.pdf",
        w0inv_plot="assets/sp4/plots/w0inv_beta.pdf",
        latent_heat_plot="assets/sp4/plots/latent_heat.pdf",

