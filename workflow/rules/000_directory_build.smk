# A rule that creates directories necessary for the pipeline

rule directory_build:
    message: "Creating directories"
    output:
        touch("directory_build.done")
    run:
        for d in directories:
            shell("mkdir -p {}".format(d))
