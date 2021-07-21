executable   = tool
arguments    = conf

universe = vanilla

output       = outfile
error        = errorfile
log          = logfile

#+RequestedAcctGroup = "group_lapp.hess"
request_cpus   = 1
request_memory = 1024
request_disk   = 1024

should_transfer_files = yes

when_to_transfer_output = ON_EXIT 

getenv = True
queue