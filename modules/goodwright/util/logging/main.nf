//
// ANSII Colours used for terminal logging
//
def log_colours(monochrome_logs) {
    Map colorcodes = [:]

    // Reset / Meta
    colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
    colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
    colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
    colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
    colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
    colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
    colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"

    // Regular Colors
    colorcodes['black']      = monochrome_logs ? '' : "\033[0;30m"
    colorcodes['red']        = monochrome_logs ? '' : "\033[0;31m"
    colorcodes['green']      = monochrome_logs ? '' : "\033[0;32m"
    colorcodes['yellow']     = monochrome_logs ? '' : "\033[0;33m"
    colorcodes['blue']       = monochrome_logs ? '' : "\033[0;34m"
    colorcodes['purple']     = monochrome_logs ? '' : "\033[0;35m"
    colorcodes['cyan']       = monochrome_logs ? '' : "\033[0;36m"
    colorcodes['white']      = monochrome_logs ? '' : "\033[0;37m"

    // Bold
    colorcodes['bblack']     = monochrome_logs ? '' : "\033[1;30m"
    colorcodes['bred']       = monochrome_logs ? '' : "\033[1;31m"
    colorcodes['bgreen']     = monochrome_logs ? '' : "\033[1;32m"
    colorcodes['byellow']    = monochrome_logs ? '' : "\033[1;33m"
    colorcodes['bblue']      = monochrome_logs ? '' : "\033[1;34m"
    colorcodes['bpurple']    = monochrome_logs ? '' : "\033[1;35m"
    colorcodes['bcyan']      = monochrome_logs ? '' : "\033[1;36m"
    colorcodes['bwhite']     = monochrome_logs ? '' : "\033[1;37m"

    // Underline
    colorcodes['ublack']     = monochrome_logs ? '' : "\033[4;30m"
    colorcodes['ured']       = monochrome_logs ? '' : "\033[4;31m"
    colorcodes['ugreen']     = monochrome_logs ? '' : "\033[4;32m"
    colorcodes['uyellow']    = monochrome_logs ? '' : "\033[4;33m"
    colorcodes['ublue']      = monochrome_logs ? '' : "\033[4;34m"
    colorcodes['upurple']    = monochrome_logs ? '' : "\033[4;35m"
    colorcodes['ucyan']      = monochrome_logs ? '' : "\033[4;36m"
    colorcodes['uwhite']     = monochrome_logs ? '' : "\033[4;37m"

    // High Intensity
    colorcodes['iblack']     = monochrome_logs ? '' : "\033[0;90m"
    colorcodes['ired']       = monochrome_logs ? '' : "\033[0;91m"
    colorcodes['igreen']     = monochrome_logs ? '' : "\033[0;92m"
    colorcodes['iyellow']    = monochrome_logs ? '' : "\033[0;93m"
    colorcodes['iblue']      = monochrome_logs ? '' : "\033[0;94m"
    colorcodes['ipurple']    = monochrome_logs ? '' : "\033[0;95m"
    colorcodes['icyan']      = monochrome_logs ? '' : "\033[0;96m"
    colorcodes['iwhite']     = monochrome_logs ? '' : "\033[0;97m"

    // Bold High Intensity
    colorcodes['biblack']    = monochrome_logs ? '' : "\033[1;90m"
    colorcodes['bired']      = monochrome_logs ? '' : "\033[1;91m"
    colorcodes['bigreen']    = monochrome_logs ? '' : "\033[1;92m"
    colorcodes['biyellow']   = monochrome_logs ? '' : "\033[1;93m"
    colorcodes['biblue']     = monochrome_logs ? '' : "\033[1;94m"
    colorcodes['bipurple']   = monochrome_logs ? '' : "\033[1;95m"
    colorcodes['bicyan']     = monochrome_logs ? '' : "\033[1;96m"
    colorcodes['biwhite']    = monochrome_logs ? '' : "\033[1;97m"

    return colorcodes
}

def dashed_line(monochrome_logs) {
    Map colors = log_colours(monochrome_logs)
    return "-${colors.dim}--------------------------------------------------------------${colors.reset}-"
}

def goodwright_logo(monochrome_logs) {
    Map colors = log_colours(monochrome_logs)
    def logo_color = colors.bblue

    String.format(
        """\n
        ${dashed_line(monochrome_logs)}
        ${logo_color}                          _               _       _     _   ${colors.reset}
        ${logo_color}     __ _  ___   ___   __| |_      ___ __(_) __ _| |__ | |_ ${colors.reset}
        ${logo_color}    / _` |/ _ \\ / _ \\ / _` \\ \\ /\\ / / '__| |/ _` | '_ \\| __|${colors.reset}
       ${logo_color}    | (_| | (_) | (_) | (_| |\\ V  V /| |  | | (_| | | | | |_ ${colors.reset}
       ${logo_color}     \\__, |\\___/ \\___/ \\__,_| \\_/\\_/ |_|  |_|\\__, |_| |_|\\__|${colors.reset}
        ${logo_color}   |___/                                   |___/           ${colors.reset}
        ${dashed_line(monochrome_logs)}
        ${colors.bgreen}${workflow.manifest.name} v${workflow.manifest.version}${colors.reset}

        ${colors.cyan}Homepage : ${workflow.manifest.homePage}${colors.reset}
        ${dashed_line(monochrome_logs)}
        """.stripIndent()
    )
}

//
// Groovy Map summarising parameters/workflow options used by the pipeline
//
def params_summary_map(workflow, params, debug) {
    // Get major nextflow workflow parameters
    def Map workflow_summary = [:]
    if (workflow.revision) {
        workflow_summary['revision'] = workflow.revision
    }
    workflow_summary['runName']      = workflow.runName
    if (workflow.containerEngine) {
        workflow_summary['containerEngine'] = workflow.containerEngine
    }
    if (workflow.container) {
        workflow_summary['container'] = workflow.container
    }
    workflow_summary['launchDir']    = workflow.launchDir
    workflow_summary['workDir']      = workflow.workDir
    workflow_summary['projectDir']   = workflow.projectDir
    workflow_summary['userName']     = workflow.userName
    workflow_summary['profile']      = workflow.profile
    workflow_summary['configFiles']  = workflow.configFiles.join(', ')

    // Get boilerplate parameters that should be present in every pipeline
    def Map params_summary = [:]
    def boilerpate_params = new LinkedHashMap()
    boilerpate_params.put("outdir", params.outdir)
    boilerpate_params.put("tracedir", params.tracedir)
    boilerpate_params.put("publish_dir_mode", params.publish_dir_mode)
    boilerpate_params.put("max_memory", params.max_memory)
    boilerpate_params.put("max_cpus", params.max_cpus)
    boilerpate_params.put("max_time", params.max_time)
    params_summary.put("Core pipeline options", boilerpate_params)

    // Output all other parameters if we have a debug flag or just non-ignored
    def other_params = new LinkedHashMap()
    def default_ignore = ['debug', 'ignore_params', 'outdir', 'tracedir', 'publish_dir_mode', 'max_memory', 'max_cpus', 'max_time', 'test_data', 'goodwright_test_data']
    def params_ignore = params.ignore_params.split(',')
    Set params_key_set = params.keySet()

    params_key_set.each {
        if(!debug && !params_ignore.contains(it) && !default_ignore.contains(it)) {
            other_params.put(it, params.get(it))
        }
        else if(debug && !default_ignore.contains(it)) {
            other_params.put(it, params.get(it))
        }
    }
    params_summary.put("Pipeline options", other_params)

    return [ 'Core Nextflow options' : workflow_summary ] << params_summary
}

//
// Get maximum number of characters across all parameter names
//
def params_maxchars(params_map) {
    Integer max_chars = 0
    for (group in params_map.keySet()) {
        def group_params = params_map.get(group)  // This gets the parameters of that particular group
        for (param in group_params.keySet()) {
            if (param.size() > max_chars) {
                max_chars = param.size()
            }
        }
    }
    return max_chars
}

//
// Beautify parameters for summary and return as string
//
def params_summary_log(workflow, params, monochrome_logs, debug) {
    Map colors = log_colours(monochrome_logs)
    String output  = ''
    def params_map = params_summary_map(workflow, params, debug)
    def max_chars  = params_maxchars(params_map)
    for (group in params_map.keySet()) {
        def group_params = params_map.get(group)  // This gets the parameters of that particular group
        if (group_params) {
            output += colors.bold + group + colors.reset + '\n'
            for (param in group_params.keySet()) {
                output += "  " + colors.blue + param.padRight(max_chars) + ": " + colors.green +  group_params.get(param) + colors.reset + '\n'
            }
            output += '\n'
        }
    }
    output += dashed_line(params.monochrome_logs)
    return output
}

//
// Generate summary log thats printable to screen
//
def summary_log(workflow, params, debug, monochrome_logs) {
    def summary_log = ''
    summary_log += goodwright_logo(monochrome_logs)
    summary_log += params_summary_log(workflow, params, monochrome_logs, debug)
    return summary_log
}

//
// Get workflow summary for MultiQC
//
def multiqc_summary(workflow, params) {
    def summary = params_summary_map(workflow, params, false)
    String summary_section = ''
    for (group in summary.keySet()) {
        def group_params = summary.get(group)  // This gets the parameters of that particular group
        if (group_params) {
            summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
            summary_section += "    <dl class=\"dl-horizontal\">\n"
            for (param in group_params.keySet()) {
                summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
            }
            summary_section += "    </dl>\n"
        }
    }

    String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
    yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
    yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
    yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
    yaml_file_text        += "plot_type: 'html'\n"
    yaml_file_text        += "data: |\n"
    yaml_file_text        += "${summary_section}"
    return yaml_file_text
}
