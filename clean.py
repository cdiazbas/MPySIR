###############################################################
#  MPySIR: MPI python script for SIR
#
#  SCRIPT: clean.py 
###############################################################

def clean():
	import os
	# We delete all node's folder
	os.system('rm -Rf node*')
	print('Directorio limpio.')
	return

if __name__ == '__main__':

	clean()
