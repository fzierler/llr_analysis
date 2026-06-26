import functools
import pandas as pd
import os

os.environ["GKSwstype"] = "100"


def list_datasets_Nt(metadata, Nt):
    return [
        f"{dataset.Nt}x{dataset.Ns}_{dataset.replicas}replicas"
        for dataset in metadata.itertuples()
        if dataset.Nt == Nt
    ]


def list_datasets(metadata):
    return [
        f"{dataset.Nt}x{dataset.Ns}_{dataset.replicas}replicas"
        for dataset in metadata.itertuples()
    ]


def list_datasets_poly_Nt(metadata,Nt):
    return [
        f"{dataset.Nt}x{dataset.Ns}_{dataset.replicas}replicas"
        for dataset in metadata.itertuples()
        if dataset.Nt == Nt
        if dataset.polyakov
    ]


def list_datasets_poly(metadata):
    return [
        f"{dataset.Nt}x{dataset.Ns}_{dataset.replicas}replicas"
        for dataset in metadata.itertuples()
        if dataset.polyakov
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
        h5file="tmp/{group}/{group}.hdf5",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file_unsorted {output.h5file} --metadata {input.metadata}'


rule sort_hdf5:
    input:
        script="scripts/sort_an.jl",
        h5file="tmp/{group}/{group}.hdf5",
        julia_instantiated="tmp/julia_ready",
    output:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file_unsorted {input.h5file} --h5file {output.h5file}'


rule tables:
    input:
        script="scripts/tables.jl",
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        julia_instantiated="tmp/julia_ready",
    output:
        table="assets/{group}/tables/runs.tex",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --outfile {output.table}'


rule overview_plots:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/trajectory_overview.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="{loc}/{group}/plots/overview/{Nt}x{Ns}_{Nreplicas}replicas.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --plot_file {output.plot} --run_name {wildcards.Nt}x{wildcards.Ns}_{wildcards.Nreplicas}replicas'


rule an_trajectory_plots:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/an_history.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="{loc}/{group}/plots/an_trajectories/{Nt}x{Ns}_{Nreplicas}replicas.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --plot_file {output.plot} --run_name {wildcards.Nt}x{wildcards.Ns}_{wildcards.Nreplicas}replicas'


rule free_energy_plots:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
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
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/compare_volumes.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/an_volume_Nt{Nt}.pdf",
    conda:
        "envs/environment.yml"
    shell:
        r'julia --project="." {input.script} --h5file {input.h5file} --Nt {wildcards.Nt} --plot_file {output.plot} --title "\$N_t={wildcards.Nt}\$"'


rule an_comparison_all:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/compare_volumes.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/an_all_runs.pdf",
    conda:
        "envs/environment.yml"
    shell:
        r'julia --project="." {input.script} --h5file {input.h5file} --plot_file {output.plot} --largets_replicas false --title ""'


rule an_replica_comparison_plots:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/compare_volumes.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/an_replicas_Nt{Nt}_Ns{Ns}.pdf",
    conda:
        "envs/environment.yml"
    shell:
        r'julia --project="." {input.script} --largets_replicas false --h5file {input.h5file} --plot_file {output.plot} --Nt {wildcards.Nt} --Ns {wildcards.Ns} --title "\$N_t\\times N_s^3={wildcards.Nt}\\times{wildcards.Ns}^3\$"'


rule entropy_plots:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/entropy.jl",
        julia_instantiated="tmp/julia_ready",
        entropy="metadata/critical_entropy_{group}.csv",
    output:
        plot="assets/{group}/plots/entropy_Nt{Nt}.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --critical_entropy {input.entropy} --plot_file {output.plot} --Nt {wildcards.Nt}'


rule critical_cumulant:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/critical_cumulants.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        csv="data_assets/{group}/critical_beta_cumulants.csv",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --outfile {output.csv}'


rule supercooling:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/supercooling.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        csv="data_assets/{group}/supercooling_param.csv",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --outfile {output.csv}'


rule critical_beta_ratio:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/critical_beta.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        csv="data_assets/{group}/critical_beta_{peak1}:{peak2}.csv",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --outfile {output.csv} --peak1 {wildcards.peak1} --peak2 {wildcards.peak2}'


rule double_gaussian_volume_plots:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/double_gaussian_volumes.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/plaquette_distribution_Nt{Nt}_volumes.pdf",
    conda:
        "envs/environment.yml"
    shell:
        r'julia --project="." {input.script} --h5file {input.h5file} --plotfile {output.plot} --Nt {wildcards.Nt} --title "\$N_t={wildcards.Nt}\$"'


rule critical_beta_volume_plots:
    input:
        script="scripts/plot_beta.jl",
        critical_beta_one_to_one="data_assets/{group}/critical_beta_1:1.csv",
        critical_beta_two_to_one="data_assets/{group}/critical_beta_2:1.csv",
        critical_beta_one_to_two="data_assets/{group}/critical_beta_1:2.csv",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/critical_beta_ratios_Nt{Nt}.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --plotfile {output.plot} --Nt {wildcards.Nt} {input.critical_beta_one_to_one} {input.critical_beta_two_to_one} {input.critical_beta_one_to_two}'


rule double_gaussian_plots_ratios:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/double_gaussian_fit.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/plaquette_distribution/pd_{Nt}x{Ns}_{Nreplicas}replicas_{peak1}:{peak2}.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --plot_file {output.plot} --peak1 {wildcards.peak1} --peak2 {wildcards.peak2} --run_name {wildcards.Nt}x{wildcards.Ns}_{wildcards.Nreplicas}replicas'


rule cumulant_plots:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/plot_cumulants.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot_binder_cumulant="assets/{group}/plots/binder_cumulant_Nt{Nt}.pdf",
        plot_specific_heat="assets/{group}/plots/specific_heat_Nt{Nt}.pdf",
    conda:
        "envs/environment.yml"
    threads: workflow.cores / 2
    shell:
        'julia --threads {threads} --project="." {input.script} --h5file {input.h5file} --plot_file_binder_cumulant {output.plot_binder_cumulant} --plot_file_specific_heat {output.plot_specific_heat} --Nt {wildcards.Nt}'


rule surface_tension_plot:
    input:
        script="scripts/surface_tension_term.jl",
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/surface_tension_term.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --plotfile {output.plot} --h5file {input.h5file}'


rule critical_beta_plot:
    input:
        script="scripts/plot_critical_beta.jl",
        csv_cumulant="data_assets/{group}/critical_beta_cumulants.csv",
        csv_histogram="data_assets/{group}/critical_beta_1:1.csv",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/beta_critical_Nt{Nt}.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --plot_file {output.plot} --input_cumulants {input.csv_cumulant} --input_histogram {input.csv_histogram} --Nt {wildcards.Nt}'


rule critical_beta_table:
    input:
        script="scripts/tex_critical_beta.jl",
        csv_cumulant="data_assets/{group}/critical_beta_cumulants.csv",
        csv_histogram="data_assets/{group}/critical_beta_1:1.csv",
        julia_instantiated="tmp/julia_ready",
    output:
        textable="assets/{group}/tables/beta_critical.tex",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --tex_file {output.textable} --input_cumulants {input.csv_cumulant} --input_histogram {input.csv_histogram} '


rule ployakov_loop_overview:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/polyakov_loop_hist.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot="assets/{group}/plots/polyakov_loop_hist/polyakov_{Nt}x{Ns}_{Nreplicas}replicas.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file {input.h5file} --plot_name {output.plot} --run_name {wildcards.Nt}x{wildcards.Ns}_{wildcards.Nreplicas}replicas'


rule ployakov_loop_vs_beta:
    input:
        h5file="data_assets/{group}/all_{group}_sorted.hdf5",
        script="scripts/polyakov_loop_moments.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        h5file="data_assets/{group}/polyakov/polyakov_loop_data.hdf5",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --h5file_in {input.h5file} --h5file_out {output.h5file}'


rule energy_histogram:
    input:
        h5file="tmp/{group}/{group}.hdf5",
        script="scripts/energy_histogram.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot_traj="assets/{group}/plots/energies/energy_trajectory_{dataset}.pdf",
        plot_hist="assets/{group}/plots/energies/energy_histogram_{dataset}.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --file {input.h5file} --name {wildcards.dataset} --plot_dir assets/{wildcards.group}/plots/energies/'


rule energy_histogram_therm:
    input:
        h5file="tmp/{group}/{group}.hdf5",
        script="scripts/energy_histogram.jl",
        julia_instantiated="tmp/julia_ready",
    output:
        plot_traj="assets/{group}/plots/energies_extra_therm/energy_trajectory_{dataset}.pdf",
        plot_hist="assets/{group}/plots/energies_extra_therm/energy_histogram_{dataset}.pdf",
    conda:
        "envs/environment.yml"
    shell:
        'julia --project="." {input.script} --file {input.h5file} --extra_therm 300 --name {wildcards.dataset} --plot_dir assets/{wildcards.group}/plots/energies_extra_therm/'


#rule ployakov_loop_susceptibility:
#    input:
#        script="scripts/polyakov_susceptibility.jl",
#        julia_instantiated="tmp/julia_ready",
#        h5files="data_assets/{group}/polyakov/polyakov_loop_data_{Nt}x{data}.hdf5"
#    output:
#        plot="assets/{group}/plots/polyakov_susceptibility_{Nt}.pdf",
#    conda:
#        "envs/environment.yml"
#    shell:
#        'julia --project="." {input.script} --plotfile {output.plot} {input.h5files}'

