###############################################################
#  MPySIR: MPI python script for SIR
#
#  SCRIPT: clean.py 
###############################################################

def clean():
	import os
	# We delete all node's folder
	os.system('rm -Rf node*')
	print('[INFO] Cleaned')
	return

if __name__ == '__main__':

	clean()
