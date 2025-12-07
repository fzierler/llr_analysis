import functools
import pandas as pd


def list_hdf5_datasets_Nt(metadata, Nt):
    return [
        f"{dataset.Nt}x{dataset.Ns}_{dataset.replicas}replicas"
        for dataset in metadata.itertuples()
        if dataset.Nt == Nt
    ]


def list_hdf5_datasets(metadata):
    return [
        f"{dataset.Nt}x{dataset.Ns}_{dataset.replicas}replicas"
        for dataset in metadata.itertuples()
    ]


@functools.cache
def parse_skip(skip):
    if skip == "[]":
        return []
    return list(map(int, skip.strip("[]").split(",")))


rule julia_instantiate:
    input:
        script="scripts/instantiate.jl",
    output:
        julia_instantiated="tmp/julia_ready",
    conda:
        "envs/environment.yml"
    shell:
        "julia {input.script} && touch {output.julia_instantiated}"


rule parse_hdf5:
    input:
        script="scripts/parse_llr.jl",
        julia_instantiated="tmp/julia_ready",
        metadata="metadata/runs_{group}.csv",
    output:
        h5file="tmp/{group}/{group}_Nt{Nt}.hdf5",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file_unsorted {output.h5file} --metadata {input.metadata} --Nt {wildcards.Nt}'


rule sort_hdf5:
    input:
        script="scripts/sort_an.jl",
        h5file="tmp/{group}/{group}_Nt{Nt}.hdf5",
        julia_instantiated="tmp/julia_ready",
    output:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file_unsorted {input.h5file} --h5file {output.h5file}'


rule tables:
    input:
        script="scripts/tables.jl",
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        julia_instantiated="tmp/julia_ready",
    output:
        table="tmp/{group}/tables/runs_Nt{Nt}.tex",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --outfile {output.table}'


rule overview_plots:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/trajectory_overview.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="tmp/{group}/plots/overview/{Nt}x{Ns}_{Nreplicas}replicas.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --plot_file {output.plot} --run_name {wildcards.Nt}x{wildcards.Ns}_{wildcards.Nreplicas}replicas'


rule an_trajectory_plots:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/an_history.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="tmp/{group}/plots/an_trajectories/{Nt}x{Ns}_{Nreplicas}replicas.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --plot_file {output.plot} --run_name {wildcards.Nt}x{wildcards.Ns}_{wildcards.Nreplicas}replicas'


rule free_energy_plots:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/free_energy.jl",
        julia_instantiated="tmp/julia_ready",
        entropy="metadata/critical_entropy_{group}.csv",
    output:
        plot="assets/{group}/plots/free_energy/free_energy_{Nt}x{Ns}_{Nreplicas}replicas.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --critical_entropy {input.entropy} --plot_file {output.plot} --run_name {wildcards.Nt}x{wildcards.Ns}_{wildcards.Nreplicas}replicas'


rule an_volume_comparison_plots:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/compare_volumes.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/an_volume_Nt{Nt}.pdf",
    conda:
        "envs/environment.yml"
    shell:
        r'julia --project="." {input.script} --h5file {input.h5file} --plot_file {output.plot} --title "\$N_t={wildcards.Nt}\$"'


rule an_replica_comparison_plots:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/compare_replicas.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/an_replicas_Nt{Nt}_Ns{Ns}.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --plot_file {output.plot} --Nt {wildcards.Nt} --Ns {wildcards.Ns}'


rule entropy_plots:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/entropy.jl",
        julia_instantiated="tmp/julia_ready",
        entropy="metadata/critical_entropy_{group}.csv",
    output:
        plot="assets/{group}/plots/entropy_Nt{Nt}.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --critical_entropy {input.entropy} --plot_file {output.plot}'


rule critical_beta:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/critical_beta.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        csv="tmp/{group}/critical_beta_Nt{Nt}.csv",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --outfile {output.csv}'


rule critical_cumulant:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/critical_cumulants.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        csv="tmp/{group}/critical_beta_cumulants_Nt{Nt}.csv",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --outfile {output.csv}'


rule critical_beta_two_to_one:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/critical_beta.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        csv="tmp/{group}/critical_beta_2:1_Nt{Nt}.csv",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --outfile {output.csv} --peak1 2'


rule critical_beta_one_to_two:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/critical_beta.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        csv="tmp/{group}/critical_beta_1:2_Nt{Nt}.csv",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --outfile {output.csv} --peak1 1 --peak2 2'


rule double_gaussian_volume_plots:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/double_gaussian_volumes.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/plaquette_distribution_Nt{Nt}_volumes.pdf",
    conda:
        "envs/environment.yml"
    shell:
        r'julia --project="." {input.script} --h5file {input.h5file} --plotfile {output.plot} --title "\$N_t={wildcards.Nt}\$"'


rule critical_beta_volume_plots:
    input:
        script="scripts/plot_beta.jl",
        critical_beta="tmp/{group}/critical_beta_Nt{Nt}.csv",
        critical_beta_two_to_one="tmp/{group}/critical_beta_2:1_Nt{Nt}.csv",
        critical_beta_one_to_two="tmp/{group}/critical_beta_1:2_Nt{Nt}.csv",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/critical_beta_volumes_Nt{Nt}.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --plotfile {output.plot} {input.critical_beta} {input.critical_beta_two_to_one} {input.critical_beta_one_to_two}'


rule double_gaussian_plots:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/double_gaussian_fit.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/plaquette_distribution/pd_{Nt}x{Ns}_{Nreplicas}replicas.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --plot_file {output.plot} --run_name {wildcards.Nt}x{wildcards.Ns}_{wildcards.Nreplicas}replicas'


rule double_gaussian_plots_two_to_one:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/double_gaussian_fit.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/plaquette_distribution/pd_{Nt}x{Ns}_{Nreplicas}replicas_two_to_one.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --plot_file {output.plot} --peak1 2 --run_name {wildcards.Nt}x{wildcards.Ns}_{wildcards.Nreplicas}replicas'


rule double_gaussian_plots_one_to_two:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/double_gaussian_fit.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/plaquette_distribution/pd_{Nt}x{Ns}_{Nreplicas}replicas_one_to_two.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --plot_file {output.plot} --peak2 2 --run_name {wildcards.Nt}x{wildcards.Ns}_{wildcards.Nreplicas}replicas'


rule cumulant_plots:
    input:
        h5file="data_assets/{group}/{group}_Nt{Nt}_sorted.hdf5",
        script="scripts/plot_cumulants.jl",
        julia_instantiated="tmp/julia_ready",
        csv="tmp/{group}/critical_beta_cumulants_Nt{Nt}.csv",
    output:
        plot_binder_cumulant="assets/{group}/plots/binder_cumulant_Nt{Nt}.pdf",
        plot_specific_heat="assets/{group}/plots/specific_heat_Nt{Nt}.pdf",
    conda:
        "envs/environment.yml"
    threads: workflow.cores / 2
    shell:
        'julia --threads {threads} --project="." {input.script} --h5file {input.h5file} --critical_values {input.csv} --plot_file_binder_cumulant {output.plot_binder_cumulant} --plot_file_specific_heat {output.plot_specific_heat} --Nt {wildcards.Nt}'


rule surface_tension_plot:
    input:
        script="scripts/surface_tension_term.jl",
        h5file_Nt4="data_assets/{group}/{group}_Nt4_sorted.hdf5",
        h5file_Nt5="data_assets/{group}/{group}_Nt5_sorted.hdf5",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/surface_tension_term.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --plotfile {output.plot} {input.h5file_Nt4} {input.h5file_Nt5}'


rule critical_beta_plot:
    input:
        script="scripts/plot_critical_beta.jl",
        csv_cumulant="tmp/{group}/critical_beta_cumulants_Nt{Nt}.csv",
        csv_histogram="tmp/{group}/critical_beta_Nt{Nt}.csv",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/beta_critical_Nt{Nt}.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --plot_file {output.plot} --input_cumulants {input.csv_cumulant} --input_histogram {input.csv_histogram} '


rule critical_beta_table:
    input:
        script="scripts/tex_critical_beta.jl",
        csv_cumulant="tmp/{group}/critical_beta_cumulants_Nt{Nt}.csv",
        csv_histogram="tmp/{group}/critical_beta_Nt{Nt}.csv",
        julia_instantiated="tmp/julia_ready",
    output:
        textable="tmp/{group}/tables/beta_critical_Nt{Nt}.tex",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --tex_file {output.textable} --input_cumulants {input.csv_cumulant} --input_histogram {input.csv_histogram} '


rule definitions:
    input:
        script="scripts/definitions.jl",
        h5file_Nt4="tmp/{group}/{group}_Nt4.hdf5",
        h5file_Nt5="tmp/{group}/{group}_Nt5.hdf5",
        julia_instantiated="tmp/julia_ready",
    output:
        definitions="assets/{group}/definitions/definitions.tex",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --outfile {output.definitions} --h5file_Nt4 {input.h5file_Nt4} --h5file_Nt5 {input.h5file_Nt5} '


rule combine_ensemble_tables:
    input:
        script="scripts/combine_tables.jl",
        file1="tmp/{group}/tables/runs_Nt5.tex",
        file2="tmp/{group}/tables/runs_Nt4.tex",
        julia_instantiated="tmp/julia_ready",
    output:
        file_out="assets/{group}/tables/runs.tex",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --outfile {output.file_out} --file1 {input.file1} --file2 {input.file2}'


rule combine_csv:
    input:
        script="scripts/combine_csv.jl",
        julia_instantiated="tmp/julia_ready",
        files_cumulants=expand(
            "tmp/{{group}}/critical_beta_cumulants_Nt{Nt}.csv", Nt=[4, 5]
        ),
        files_critical=expand("tmp/{{group}}/critical_beta_Nt{Nt}.csv", Nt=[4, 5]),
        files_ratios=expand(
            "tmp/{{group}}/critical_beta_{r}_Nt{Nt}.csv", Nt=[4, 5], r=["1:2", "2:1"]
        ),
    output:
        out_cumulants="data_assets/{group}/critical_beta_cumulants.csv",
        out_critical="data_assets/{group}/critical_beta_distribution.csv",
        out_ratios="data_assets/{group}/critical_beta_ratios.csv",
    conda:
        "envs/environment.yml"
    shell:
        """
        julia --project="." {input.script} --outfile {output.out_cumulants} {input.files_cumulants}
        julia --project="." {input.script} --outfile {output.out_critical} {input.files_critical}
        julia --project="." {input.script} --outfile {output.out_ratios} {input.files_ratios}
        """


rule combine_critical_beta_tables:
    input:
        script="scripts/combine_tables.jl",
        file1="tmp/{group}/tables/beta_critical_Nt5.tex",
        file2="tmp/{group}/tables/beta_critical_Nt4.tex",
        julia_instantiated="tmp/julia_ready",
    output:
        file_out="assets/{group}/tables/beta_critical.tex",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --outfile {output.file_out} --file1 {input.file1} --file2 {input.file2}'
