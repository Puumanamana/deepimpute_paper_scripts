data_dir = "${workflow.projectDir}/../paper_data/speed_and_memory"

process RscriptRunner {
    tag { "${method}_${ncells}" }
    publishDir "Results/${method}", mode: "copy"
    
    input:
    each ncells from params.ncell_list
    each trial from params.trials
    each method from params.methods_R
    
    output:
    file "{vmstats, time}*"
    
    script:
    """
    nohup vmstat 1 10000 > vmstats_${method}_${ncells}_${trial} &
    Rscript ${workflow.projectDir}/imputation_runner.R ${data_dir}/mouse1M_${ncells}_nonTransposed.csv ${method} ${trial} ${task.cpus}
    pkill vmstat
    """
}

process PythonRunner {
    tag { "${method}_${ncells}" }
    publishDir "Results/${method}", mode: "copy"
    
    input:
    each ncells from params.ncell_list
    each trial from params.trials
    each method from params.methods_py
    
    output:
    file "{vmstats, time}*"
    
    script:
    """
    nohup vmstat 1 10000 > vmstats_${method}_${ncells}_${trial} &
    python3 ${workflow.projectDir}/imputation_runner.py \
         --path ${data_dir}/mouse1M_${ncells}_transposed.csv \
         --method ${method} \
         --threads ${task.cpus} \
         --trial ${trial}
    pkill vmstat
    """
}
