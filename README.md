## NGS-MSTB: Next Generation Sequencing Microbial Surveillance Toolbox for scalable de-novo assembly of viral genomes and bacterial genes

This repository MICGENT is the main component of the NGS-MSTB system.
Here we describe how to use the entire system that is comprised of
multiple code repositories.

## Description

### Summary 

NGS-MSTB is a bioinformatics pipeline for quickly and reliably assembling next generation sequencing reads from large scale viral and bacterial clinical strain surveillance studies to monitor potential drug resistance and the presence of virulence factors. Using the pipeline’s integrated graphical interface, an analyst can routinely assemble and review thousands of datasets in one batch within a few hours on a high-performance compute cluster.

The software assembles de-novo complete viral genomes from amplicon-derived shotgun reads potentially sequenced with extremely high local coverage variability as well as specific bacterial genes from whole genome shotgun reads.

The software was successfully deployed as the assembly module in a system validated according to good clinical practices (GCP). To this end, it guarantees deterministic results, full provenance data and ensures the integrity of the inputs and intermediate datasets through cryptographic checksums.

### Availability and Implementation

Source code, binary Conda packages and testing datasets are freely available from this and several associated repositories on GitHub as well as on Zenodo. The pipeline provides a user interface based on the Galaxy bioinformatics workbench, dynamic Web reports with genome browser views and a distributed backend coded with Common Workflow Language and Python. A Docker container is available to simplify deployment.

## Author

Andrey Tovchigrechko `<andreyto AT gmail.com>`

## License

GPLv3. See also COPYING file that accompanies the source code.

## Citation

To learn about the motivation and the algorithm please see the manuscript:
bioRxiv URL

If you find this software useful, please cite the manuscript.

## How to run - the overview

Unless you need a very fast turnaround (a few hours) on batches of thousands of sequencing datasets (samples), you will probably be OK with running our provided [NGS-MSTB Docker image](https://hub.docker.com/orgs/ngsmstb) on a single multicore host machine. The machine can be a virtual machine (VM) in a commercial cloud, on-premise server or a laptop. Increasing the number of cores and, proportionally, RAM on the host machine will give you a roughly linear increase in performance over the number of samples. This is generally known as "vertical scaling". Using a Docker image greatly simplifies deployment of this complex integrated pipeline. For scaling beyond a single multicore machine, you will have to deploy the NGS-MSTB Conda packages and Galaxy tools into the existing Galaxy instance that is configured to run jobs on a distributed compute cluster comprising multiple nodes (horizontal scaling).

Currently, the dynamic Web report that NGS-MSTB generates will load several open-source JavaScript libraries from the public Content Distribution Networks (CDNs). This means that the client machine where you view the report in a Web browser will need access to the Internet.

## How to run - Docker

This assumes that you already have got some system operational on the host machine allowing you to run Docker containers. Our examples use Docker commands entered from a Linux terminal. You can easily adapt the commands to a different container execution system such as Redhat Podman and/or different host operating system (MacOS or Windows). Most of our testing of the NGS-MSTB container has been done
on a Linux VM with 64GB RAM and 16 CPU cores and Mac PowerBook Pro with 32GB RAM
and 6 CPU cores (12 cores with hyperthreading)

### The minimum hardware resources

The Docker containers should be able to use at least four CPU cores and 10 GB of RAM. On Linux hosts, the containers by default have access to all CPUs and RAM of the host. On Mac and Windows, the containers run inside a VM, and the default Docker settings restrict the VM resource allocation. You should check and adjust these if necessary. For example, under Docker Desktop, these settings will be found under _Preferences -> Advanced_.
See disk space requirements below.

If you are increasing the hardware resources in order to process your data faster, make sure that you **increment cores and RAM proportionally**: for each four additional CPU cores, add eight GB of RAM.

### Quick Start

### Start the container:

`docker run --rm -d -p 8080:80 --name ngs-mstb -e GALAXY_CONFIG_SINGLE_USER=admin_ge@ngs-mstb.nowhere ngsmstb/ngs-mstb:latest`

### Interact with Galaxy:

In your Web browser, open the URL `localhost:8080`. If you started the container on a remote machine, replace `localhost` in the URL with the DNS name of that machine. You can also use a different port instead of `8080` both in the `docker run` command and in the URL. The recommended Web browser is Chrome, Firefox or other standards-complying browser.

The Galaxy instance inside the container will take about a minute to become fully operational - you might see `Connection refused` when you first type the URL and have to refresh the page a couple of times.

Eventually, you will see the Galaxy user interface (UI) in your browser window. Follow the links at the center to take the interactive tours of executing the example of assembling several complete genomes of Respiratory Syncytial Virus (RSV) and reviewing the resulting Web report. All necessary example datasets are included in the container.

### Required disk space

The Quick Start command above will use transient anonymous Docker volumes for both temporary working files during pipeline execution and the final outputs stored as Galaxy datasets. This has several implications:
- Once you stop the container (`docker stop`) or reboot the host, 
  your output data and all Galaxy histories will disappear
- The file system where your Docker installation creates volumes should 
  have enough free space to hold the temporary pipeline data. The pipeline caches all intermediate output files during its run for 
  the ability to automatically retry  any failed workflow steps. This means that it will use a total of 10 GB at peak time in the case of the eight RSV genomes example running on four CPU cores. Higher degree of parallelism (more CPU cores in the container) will increase peak disk use because more transient work
  files will coexist in the workflow at any given time. On Mac and Windows, Docker volumes typically reside on the same file 
  system as the host operating system (OS). On Linux, the volumes might be on a relatively small separate partition that can 
  overfill easily.

For the use beyond running the Quick Start example, you might want to supply a "bind-mount" that will make your outputs persist 
across container restarts and allow using any location you want for the pipeline temporary data:

```
docker run --rm -d -p 8080:80 --name ngs-mstb \
-v "/host/directory/path":"/export" \
-e GALAXY_CONFIG_SINGLE_USER=admin_ge@ngs-mstb.nowhere ngsmstb/ngs-mstb:latest
```

Above, we used the Unix shell line continuation symbol `\` to split the long
command across three lines. You should remove this symbol if concatenating
the command into a single line.

In the command above, the placeholder `/host/directory/path` should be replaced with the actual absolute path to the existing directory on your 
host machine. That directory should be empty when you first run this command and be located on a file system with enough free space. Galaxy will copy multiple files into that directory when it first starts. If you restart the container and use the same command again, Galaxy will notice the existing files and restore your previous session.

**Note**: Your Docker installation will still need to have enough free space to keep the container image itself (about 6GB). If you try a lot of different containers, this space fills up quickly, which creates a common source of problems encountered while using Docker in general. The housekeeping commands
such as `docker system prune` and `docker volume prune` in case of Docker 
might help in such situations.

### How to supply your NGS sequencing inputs

MGS-MSTB needs two types of externally provided data files:

- FASTQ files with paired-end NGS reads provided as two separate files for
  forward and reverse reads per each sequencing sample.
- One multi-FASTA file with recruiting references. Depending on the 
  the analysis objectives, the same file can be used as QC references for
  aligning assembly contigs, or a second multi-FASTA file can be supplied
  for this purpose.

#### FASTQ read files

The FASTQ files should be named in a consistent pattern to allow configuring
a regular expression for extracting sample IDs. These IDs will be used throughout
the pipeline to name rows in various output tables and to name the assembly
contigs in the output FASTA files.

Let us assume that you have placed the FASTQ files under some directory named 
`/path/to/seqstore/data` on the host file system where you start your Docker container.

Then, you should bind-mount this directory under a specific location 
`/seqstore/data` in the container:

```
docker run --rm -d -p 8080:80 --name ngs-mstb \
-v "/path/to/seqstore/data":"/seqstore/data" \
-v "/host/directory/path":"/export" \
-e GALAXY_CONFIG_SINGLE_USER=admin_ge@ngs-mstb.nowhere ngsmstb/ngs-mstb:latest
```

The container will only need read-only access to that location. **Note**:
the pipeline in the container will be executing as user `galaxy` that
has user ID (UID) 1450 and the same group ID (GID). You should either
create a user/group with the same UID/GID and give read and directory
execute permissions to either user or group to all files under your 
`/path/to/seqstore/data`, or make `/path/to/seqstore/data` readable by all
users (if that is admissible for you security-wise). The latter can be
done with the following command on Linux or MacOS:
`chmod -R o+rX "/path/to/seqstore/data"`. Replace `/path/to/seqstore/data` everywhere in these instructions with your actual absolute path.

Pairs of FASTQ files across multiple samples can be spread across several subdirectories
under `/path/to/seqstore/data` - please see the Help lines in the 
manifest builder tool in NGS-MSTB Galaxy container. The Help describes both
multi-directory search patterns and the regular expressions for extracting
the sample IDs.

When specifying path to the data subdirectory in the NGS-MSTB manifest building
tool, the **path should be given as relative to the fixed root path `/seqstore` inside the container**.

This is easier to demonstrate through the example:

Let us assume that on your host machine the FASTQ files are located under
a subdirectory `/path/to/seqstore/data/my_sequencing_run1`, and their names end
in `.fastq.gz` Then, in the manifest tool, you would supply:
`data/my_sequencing_run1/*.fastq.gz`.

#### FASTA reference files

These should be uploaded into your Galaxy history through the Web
browser using Galaxy Get Data tool. Please see the interactive tour
in the NGS-MSTB Galaxy container.

### Multi-user Galaxy instance

In the commands above, the environment variable `-e GALAXY_CONFIG_SINGLE_USER=admin_ge@ngs-mstb.nowhere` instructs Galaxy to starts in a single-user mode and automatically login the default Galaxy administrator account. Although you can use Galaxy UI in an anonymous session without logging in, you will not get any persistence of your work history if your restart the container or even close the browser window. This is why we used the automatic login above. However, this is only acceptable if your Galaxy container is executed on an isolated host and not shared between several user. If you want to use it in a shared environment, you should at least protect the Galaxy admin access by overriding the default values 
of several access keys:

```
docker run --rm -d -p 8080:80 --name ngs-mstb -v "/host/directory/path":"/export" \
-e GALAXY_CONFIG_MASTER_API_KEY=your_unique_master_key \
-e GALAXY_DEFAULT_ADMIN_PASSWORD=your_unique_admin_password \
-e GALAXY_DEFAULT_ADMIN_KEY=your_unique_admin_key \
-e NGS_MSTB_SIG_KEY=your_unique_signature_key \
-e NGS_MSTB_GALAXY_REPORTS_PASSWORD=your_unique_galaxy_reports_password \
ngsmstb/ngs-mstb:latest
```

Replace each `your_unique_*` placeholder above with some secret string. The default
values of all those variables baked into the Docker image can be found inside the 
`Dockerfile` in the `ngs-mstb-docker` repository or printed by inspecting the image
with `docker inspect ngsmstb/ngs-mstb:latest`. In particular, default admin user 
is `admin_ge@ngs-mstb.nowhere` with a password `NgsMstb20`.

You can login as the admin user and access the standard Galaxy `Admin` menu to create additional users for a shared deployment scenario.

You can look at the section "How to run - Docker for a validated GCP deployment" for more complete set of options of configuring production secure deployments of NGS-MSTB.

## How to run - Docker for a validated GCP deployment

As part of `ngs-mstb-docker` repository, we have provided the deployment 
scripts which pull the key configuration parameters from another repository
called `ngs-mstb-secrets`. The latter is just the example. The overall design
assumes that you would create a private copy of the `ngs-mstb-secrets`, tag it and
populate it with your own secret values for the critical access keys, Galaxy 
Web page branding, certificates for HTTPs encrypted communication and
LDAP authentication config. These can be customized inside the same repository to configure several target hosts
on which you would want to deploy the NGS-MSTB container. That allows easily
deploying on the development, testing and production VMs as it is typical in
a validated software environment.

Under the `ngs-mstb-docker` repository, we have provided an example
set of instructions in a file `example_prepare_host.sh`. These instructions demonstrate how to prepare the Linux host
running Red Hat Enterprise Linux for deploying the NGS-MSTB container.
The instructions include creating the matching accounts in the host OS
and setting proper file permissions on the data directories.

Once the host was prepared, the script `docker_deploy.sh` can be executed.
This script will pull the customized configs and other run-time artifacts
from the `ngs-mstb-secrets` for the target host and launch the container
using the Docker Compose service definition. The latter will ensure that
the container will get automatically restarted with the same configuration
every time the host is rebooted.

## In depth understanding of the NGS-MSTB Galaxy tools

The NGS-MSTB Galaxy tools provide extensive Help strings linked to all parameters, as well as 
tool-level Help sections located below the Execute button on each tool form.
Those tool-level Help sections describe what the tool does, its place in the
overall workflow and the outputs.

## More options for running the container

The NGS-MSTB Docker image is based on the [Galaxy Docker image](https://github.com/bgruening/docker-galaxy-stable). The upstream Web page describes
in detail various additional configuration options which can be used at
run-time.

## Automated testing suite

The software includes a testing suite that is composed both from unit tests and integration tests.
The integration tests include checks for determinism. Those are implemented by running the same
samples in multiple replicates and in several batches, and then comparing outputs for complete
identity across replicates and batches.

The anonymized input sequencing reads used in the tests are available from 
[Zenodo NGS-MSTB community page](https://zenodo.org/communities/ngs-mstb).
In the repository `ngs-mstb-docker`, the script `docker_run_tests.sh` will execute the tests
inside a container.
