include: "rules.smk"


metadata = pd.read_csv("metadata/runs_su3_zeroT.csv")


rule all:
    input:
        overview=expand(
            "assets/su3_zeroT/plots/overview/{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        trajectory=expand(
            "assets/su3_zeroT/plots/an_trajectories/{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        an_volume=expand("assets/su3_zeroT/plots/an_volume_Nt{Nt}.pdf", Nt=[4, 12, 16]),
        an_all="assets/su3_zeroT/plots/an_all_runs.pdf",
        an_replicas="assets/su3_zeroT/plots/an_replicas_Nt12_Ns12.pdf",
