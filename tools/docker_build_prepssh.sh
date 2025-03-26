mkdir -p ~/.ssh
touch ~/.ssh/ldg_user ~/.ssh/ldg_token
echo "${OSG_ACCESS}" > ~/.ssh/id_rsa
echo -e "Host sugwg-test1.phy.syr.edu\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config ;
echo -e "Host sugwg-condor.phy.syr.edu\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config ;
echo -e "Host oasis-login.opensciencegrid.org\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config ;
echo -e "Host code.pycbc.phy.syr.edu\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config ;
echo -e "@cert-authority * ecdsa-sha2-nistp256 AAAAE2VjZHNhLXNoYTItbmlzdHAyNTYAAAAIbmlzdHAyNTYAAABBBHa03AZF3CvJ1C4Po15swSaMYI4kPszyBH/uOKHQYvu+EpehSfMZMaX5D7pUpc5cAXvMEEFzlZJQH4pOioIlqyE= IGWN_CIT_SSH_CERT_AUTHORITY" >> ~/.ssh/known_hosts
chmod 600 ~/.ssh/id_rsa ~/.ssh/config ~/.ssh/ldg_user ~/.ssh/ldg_token
