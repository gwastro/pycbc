mkdir -p ~/.ssh
touch ~/.ssh/id_rsa ~/.ssh/ldg_user ~/.ssh/ldg_token
chmod 600 ~/.ssh/id_rsa ~/.ssh/ldg_user ~/.ssh/ldg_token
echo "${OSG_ACCESS}" > ~/.ssh/id_rsa
echo -e "Host sugwg-test1.phy.syr.edu\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config ;
echo -e "Host sugwg-condor.phy.syr.edu\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config ;
echo -e "Host oasis-login.opensciencegrid.org\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config ;
echo -e "Host code.pycbc.phy.syr.edu\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config ;
chmod 600 ~/.ssh/id_rsa ~/.ssh/config ~/.ssh/ldg_user ~/.ssh/ldg_token
