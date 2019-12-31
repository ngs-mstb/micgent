To prepare package list for major Conda upgrade:
`source activate ariba`
`conda list --export > ariba-env.txt`
Then, perform regex replacement to remove version numbers.
To take a snapshot of installed versions including pip packages:
`conda env export > conda_env_ariba.yml`

