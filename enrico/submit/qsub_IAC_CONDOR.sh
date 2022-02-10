executable   = tool
arguments    = conf

universe = vanilla

output       = outfile
error        = errorfile
log          = logfile

+PreCmd        = "/bin/bash" 
+PreArguments  = "/net/nas/proyectos/magic/magic/mnievas/bin/load_fermi.sh"

#+RequestedAcctGroup = "group_lapp.hess"
request_cpus   = 1
request_memory = 1024
request_disk   = 1024
concurrency_limits = mnievas:83

#transfer_input_files = /net/nas/proyectos/magic/magic/mnievas/bin/load_fermi.sh,enricodir
transfer_input_files  = /bin/bash,/net/nas/proyectos/magic/magic/mnievas/bin/load_fermi.sh
should_transfer_files = yes

when_to_transfer_output = ON_EXIT 

getenv = True
queue
