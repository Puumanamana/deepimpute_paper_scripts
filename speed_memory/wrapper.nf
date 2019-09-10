data_dir = "${workflow.projectDir}/../paper_data/speed_and_memory"
  //methods_py = ('deepImpute','magic','dca')
  // methods_R = ('SAVER','scImpute','VIPER','DrImpute')
methods_py = ['magic']
methods_R = ['SAVER']
trials = [1,2]
ncell_list=[100,500]


process RscriptRunner {
    tag { "${method}_${ncells}" }
    publishDir "Results/${method}", mode: "copy"
    
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
    tag { "${method}_${ncells}" }
    publishDir "Results/${method}", mode: "copy"
    
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
    publishDir "Results/Summary", mode: "copy"
    
    input:
    file f from DATA_R.collect()
    file g from DATA_PY.collect()
    
    output:
    file "summary.csv"
    
    script:
    """
    cat time* > summary.csv
    for f in `ls vmstat*`; do python ${workflow.projectDir}/extract_memory.py \$f ; done
    """
}
