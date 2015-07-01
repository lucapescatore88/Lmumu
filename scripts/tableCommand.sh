cat correl_signi.txt | grep variable | grep -v Keep | sed s/":"/" "/g | awk '{ print $3, "\t& "$4, "\t& "$5, "\t& "$6, "\t& "$7, " \\\\" }'
