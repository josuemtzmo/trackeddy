import socket
import ftplib
from ftplib import FTP
import sys

#if __name__ == "__main__":
#    login = str(sys.argv[1])
#    passwd = str(sys.argv[2])
#    print(login,passwd)
#ftp = FTP('ftp.debian.org')     # connect to host, default port
#ftp.login()  
#print(ftp)
ip='ftp://ftp.sltac.cls.fr/Core/SEALEVEL_GLO_PHY_L4_REP_OBSERVATIONS_008_047/dataset-duacs-rep-global-merged-allsat-phy-l4-v3'
#port=990


ftp = FTP(ip)
ftp.login(login,passwd)




