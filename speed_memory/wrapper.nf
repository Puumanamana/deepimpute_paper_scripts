data_dir = "${workflow.projectDir}/../paper_data/speed_memory"
result_dir = "${workflow.projectDir}/../results/speed_memory/nextflow"

methods_py = ['magic','deepImpute','dca']
methods_R = ['SAVER','scImpute','VIPER','DrImpute']

trials = [1,2,3,4,5]
ncell_list=[100,500,1000,5000]

process RscriptRunner {
    tag { "${method}_${ncells}_${trial}" }
    publishDir "${result_dir}/${method}", mode: "copy"
    
    input:
    each ncells from ncell_list
    each trial from trials
    each method from methods_R
    
    output:
    file "{vmstats,time}*" into DATA_R
    
    script:
    """
    nohup vmstat 1 10000 > vmstats_${method}_${ncells}_${trial}.txt &
    Rscript ${workflow.projectDir}/imputation_runner.R ${data_dir}/mouse1M_${ncells}_nonTransposed.csv ${method} ${ncells} ${trial} ${task.cpus}
    pkill vmstat
    """
}

process PythonRunner {
    tag { "${method}_${ncells}_${trial}" }
    publishDir "${result_dir}/${method}", mode: "copy"
    
    input:
    each ncells from ncell_list
    each trial from trials
    each method from methods_py
    
    output:
    file "{vmstats,time}*" into DATA_PY
    
    script:
    """
    nohup vmstat 1 10000 > vmstats_${method}_${ncells}_${trial}.txt &
    python3 ${workflow.projectDir}/imputation_runner.py \
         --path ${data_dir}/mouse1M_${ncells}_transposed.csv \
         --method ${method} \
         --threads ${task.cpus} \
         --trial ${trial} \
         --ncells ${ncells}
    pkill vmstat
    """
}

process SummaryFiles {
    publishDir "${result_dir}/summary", mode: "copy"
    
    input:
    file f from DATA_R.collect()
    file g from DATA_PY.collect()
    
    output:
    file "time_memory.csv"
    
    script:
    """
    cat time* > time_memory.csv
    for f in `ls vmstat*`; do python ${workflow.projectDir}/extract_memory.py \$f ; done
    """
}
