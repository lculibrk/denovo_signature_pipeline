import glob
import pandas as pd

configfile: "config.yaml"

target_functions = {}
project_cohorts = {}
project_genome = {}
maf_colnames = {}

newline_re = re.compile(r'\s*\n\s*')
def replace_newlines(input_string):
    return newline_re.sub(' ', input_string)

PROJECT_ROOT = config["PROJECT_ROOT"]
WTSI_ITERATIONS_PER_CORE = config["WTSI_ITERATIONS_PER_CORE"]
WTSI_CORES = config["WTSI_CORES"]
WTSI_NSIG_RANGE = config["WTSI_NSIG_RANGE"]
REFERENCE_SIGNATURE_PATH = config["REFERENCE_SIGNATURE_PATH"]
PROJECTS = config["PROJECT_NAMES"]
PROJECT_GENOMES = config["PROJECT_GENOMES"]
MAF_COLNAMES = config["MAF_COLNAMES"]

matlab_source_paths = '''
    addpath('scripts/pipelines/wtsi-signatures');
    addpath('scripts/pipelines/wtsi-signatures/custom');
    addpath('scripts/pipelines/wtsi-signatures/source');
'''



def glob_to_df(glob_string, regex, column_names):
    regex_obj = re.compile(regex)
    glob_output = glob.iglob(glob_string)
    regex_searches = (regex_obj.search(s) for s in glob_output)
    regex_matches = [[r.group(i) for i in range(1, len(column_names)+1)] for r in regex_searches]
    df = pd.DataFrame(
        regex_matches
    )
    df = df.iloc[:, range(len(column_names))]
    df.columns = column_names
    return df


def get_snv_targets(wildcards):
    df = glob_to_df(
        'data/{project}/{method}/cohorts/{cohort}/snv/*.maf'.format(
            method = wildcards.method,
            cohort = wildcards.cohort,
	    project = wildcards.project
        ),
        r'data\/' + wildcards.project + '\/' + wildcards.method + '\/cohorts\/' + wildcards.cohort + '\/snv\/(.*?).maf',
        ['samples']
    )
    return df

def get_aggregation_targets(wildcards):
    ''' Receives the following wildcards:
            - project
            - cohort
            - mutation_type
            - method
    '''
    targets_df = get_snv_targets(wildcards)
    targets = 'output/' + \
        wildcards.project + '/' + \
        wildcards.method + '/cohorts/' + \
        wildcards.cohort + '/' + \
        wildcards.mutation_type + '/catalog/' + \
        targets_df.samples + '.tsv'
    return targets

def get_individual_aggregation_targets(wildcards):
    target_loading_function = \
        target_functions[wildcards.project][wildcards.mutation_type]
    targets_df = target_loading_function(wildcards)
    targets = 'output/' + \
        wildcards.project + '/' + \
        wildcards.method + '/cohorts/' + \
        wildcards.cohort + '/signatures/' + \
        wildcards.mutation_type + '/' + \
        wildcards.signature_method + '/' + \
        wildcards.summary_type + '/' + \
        targets_df.samples + '.tsv'
    return targets


## Create all outputs

project_outs = expand('output/{project}/genome/cohorts/{{cohort}}/signatures/snv/signature_report.html', project = PROJECTS)

#print(PROJECTS)

inputs = []
for project in PROJECTS:
    print(project)
    print(config["ANALYSIS_COHORTS"])
#    print(expand('output/{project}/genome/cohorts/{cohort}/signatures/snv/signature_report.html', cohort = config["ANALYSIS_COHORTS"][project], project = project))
    inputs.extend(expand('output/{project}/genome/cohorts/{cohort}/signatures/snv/signature_report.html', cohort = config["ANALYSIS_COHORTS"][project], project = project))

#print(inputs)
## Edit the below ALLCAPS stuff to match your data
rule all:
    input:
        inputs


########################################
### Aggregation of Mutation Catalogs ###
########################################

rule snv_catalogs:
    input:
        'data/{project}/{method}/cohorts/{cohort}/snv/{sample}.maf'
    output:
        'output/{project}/{method}/cohorts/{cohort}/snv/catalog/{sample}.tsv'
    params:
        genome = lambda wildcards: PROJECT_GENOMES[wildcards.project],
        colnames = lambda wildcards: MAF_COLNAMES[wildcards.project]
    shell:
        'Rscript scripts/pipelines/nnls-signatures/mutation_catalog.R -t {input} -o {output} --colnames {params.colnames} --genome {params.genome}'

rule catalog_paths:
    input:
        get_aggregation_targets
    output:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/catalog_paths.tsv'
    wildcard_constraints:
        summary_type='(summary|summary_fraction)'
    run:
        f = open(output[0], 'w')
        f.write('\n'.join(input))
        f.close()

rule aggregate_snv_catalogs:
    input:
        'output/{project}/{method}/cohorts/{cohort}/signatures/snv/catalog_paths.tsv'
    output:
        'output/{project}/{method}/cohorts/{cohort}/signatures/snv/catalogs.tsv'
    shell:
        'Rscript scripts/pipelines/basic-functions/row_bind_tables.R -i {input} -o {output}'

########################################
### WTSI Mutation Signature Pipeline ###
########################################

rule create_wtsi_input:
    input:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/catalogs.tsv'
    output:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/input.mat'
    shell:
        replace_newlines('''Rscript scripts/pipelines/wtsi-signatures/catalog_tidy_to_mat.R
            -i {input}
            -o {output}
            -s {wildcards.project}_{wildcards.method}_{wildcards.cohort}
            -c paths
        ''')

rule prepare_wtsi_iterations:
    input:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/input.mat'
    output:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/iteration_input.mat'
    shell:
        replace_newlines('''octave --eval "
            {matlab_source}
            outputIterationData('{{input}}', '{{output}}');
            exit
            "
        '''.format(matlab_source = matlab_source_paths))

rule run_wtsi_iterations:
    input:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/iteration_input.mat'
    output:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/cores/output_{nsig}-signatures_core-{core}.mat'
    threads: 1
    shell:
        replace_newlines('''octave --eval "
            pkg load statistics;
            {matlab_source}
            runIterations(
              '{{input}}',
              '{{output}}',
              '{iter_per_core}', '{{wildcards.nsig}}'
            );
            exit
            "
        '''.format(
            matlab_source = matlab_source_paths,
            iter_per_core = WTSI_ITERATIONS_PER_CORE
        ))

rule wtsi_iteration_paths_file:
    input:
        expand('output/{{project}}/{{method}}/cohorts/{{cohort}}/signatures/{{mutation_type}}/wtsi/cores/output_{{nsig}}-signatures_core-{core}.mat',
            core = range(WTSI_CORES)
        )
    output:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/iteration_paths/{nsig}_core_output_paths.txt'
    run:
        f = open(output[0], 'w')
        f.write('\n'.join(input))
        f.close()

rule merge_wtsi_iterations:
    input:
        original = 'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/input.mat',
        iter_paths = 'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/iteration_paths/{nsig}_core_output_paths.txt'
    output:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/output/{nsig}-signatures.mat'
    threads: 1
    shell:
        replace_newlines('''
            bash
                scripts/pipelines/wtsi-signatures/compiled/mergeIterations/for_redistribution_files_only/run_mergeIterations.sh
                /projects/ezhao_prj/software/MCR/v901
                {input.original}
                {input.iter_paths}
                {output}
        ''')

rule export_wtsi_metrics:
    input:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/output/{nsig}-signatures.mat'
    output:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/metrics/{nsig}-signatures.tsv'
    shell:
        'Rscript scripts/pipelines/wtsi-signatures/export_wtsi_metrics.R -i {input} -o {output}'

rule aggregate_wtsi_metrics:
    input:
        expand('output/{{project}}/{{method}}/cohorts/{{cohort}}/signatures/{{mutation_type}}/wtsi/metrics/{nsig}-signatures.tsv', nsig=WTSI_NSIG_RANGE)
    params:
        comma_separated = lambda wildcards, input: ','.join(input)
    output:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/model_selection_metrics.tsv'
    shell:
        'Rscript scripts/pipelines/basic-functions/row_bind_tables.R -p {params.comma_separated} -o {output}'

rule project_aggregation_signatures:
    input:
        lambda wildcards: expand(
            'output/{project}/{method}/cohorts/{{cohort}}/signatures/{mutation_type}/wtsi/combined_signatures.tsv'.format(
                project = wildcards.project,
                method = wildcards.method,
                mutation_type = wildcards.mutation_type,
            ),
            cohort=project_cohorts[wildcards.project],
        )
    params:
        comma_separated = lambda wildcards, input: ','.join(input)
    output:
        'output/{project}/{method}/aggregated/{mutation_type}/signatures.tsv'
    shell:
        'Rscript scripts/pipelines/basic-functions/row_bind_tables.R -p {params.comma_separated} -o {output}'

rule project_aggregation_exposures:
    input:
        lambda wildcards: expand(
            'output/{project}/{method}/cohorts/{{cohort}}/signatures/{mutation_type}/wtsi/combined_exposures.tsv'.format(
                project = wildcards.project,
                method = wildcards.method,
                mutation_type = wildcards.mutation_type,
            ),
            cohort=project_cohorts[wildcards.project],
        )
    params:
        comma_separated = lambda wildcards, input: ','.join(input)
    output:
        'output/{project}/{method}/aggregated/{mutation_type}/exposures.tsv'
    shell:
        'Rscript scripts/pipelines/basic-functions/row_bind_tables.R -p {params.comma_separated} -o {output}'

rule project_aggregation_metrics:
    input:
        lambda wildcards: expand(
            'output/{project}/{method}/cohorts/{{cohort}}/signatures/{mutation_type}/wtsi/model_selection_metrics.tsv'.format(
                project = wildcards.project,
                method = wildcards.method,
                mutation_type = wildcards.mutation_type,
            ),
            cohort=project_cohorts[wildcards.project],
        )
    params:
        comma_separated = lambda wildcards, input: ','.join(input)
    output:
        'output/{project}/{method}/aggregated/{mutation_type}/model_selection_metrics.tsv'
    shell:
        'Rscript scripts/pipelines/basic-functions/row_bind_tables.R -p {params.comma_separated} -o {output}'

rule export_wtsi_snv_signatures_exposures:
    input:
        expand('output/{{project}}/{{method}}/cohorts/{{cohort}}/signatures/{{mutation_type}}/wtsi/output/{nsig}-signatures.mat', nsig=WTSI_NSIG_RANGE)
    params:
        comma_separated = lambda wildcards, input: ','.join(input)
    output:
        signatures='output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/combined_signatures.tsv',
        exposures='output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/combined_exposures.tsv'
    shell:
        'Rscript scripts/pipelines/wtsi-signatures/export_wtsi_signatures_exposures.R -p {params.comma_separated} -s {output.signatures} -e {output.exposures}'

rule wtsi_snv_pipeline:
    input:
        signatures='output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/combined_signatures.tsv',
        exposures='output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/combined_exposures.tsv',
        metrics='output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/model_selection_metrics.tsv',
    output:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/COMPLETE'
    shell:
        'touch {output}'


########################
### Signature Report ###
########################

rule snv_signature_report:
    input:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/wtsi/COMPLETE'
    output:
        'output/{project}/{method}/cohorts/{cohort}/signatures/{mutation_type}/signature_report.html'
    params:
        root_dir = lambda wildcards, output: PROJECT_ROOT + '/' + '/'.join(output[0].split('/')[0:-1]),
        filename = lambda wildcards, output: output[0].split('/')[-1]
    shell:
        '''Rscript -e "rmarkdown::render('scripts/pipelines/reports/mutation_signature.Rmd', knit_root_dir='{params.root_dir}', output_dir='{params.root_dir}', output_file='{params.filename}')"'''