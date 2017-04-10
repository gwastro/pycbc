#!/bin/bash

set -e

echo -e "\\n>> [`date`] Testing pycbc_make_coinc_search_workflow"

TRAVIS_TAG=${1}

if [ "x${TRAVIS_TAG}" == "x" ] ; then
  echo -e "\\n>> [`date`] Error: TRAVIS_TAG was not passed to script or is empty"
  exit 1
fi

echo -e "\\n>> [`date`] Entering virtual environment pycbc-${TRAVIS_TAG}"
source /cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/x86_64_rhel_7/virtualenv/pycbc-${TRAVIS_TAG}/bin/activate

echo -e "\\n>> [`date`] Cloning pycbc-config git repository"
test -r pycbc-config || git clone --depth 1 git@code.pycbc.phy.syr.edu:ligo-cbc/pycbc-config.git
CONFIG_PATH="file://`pwd`/pycbc-config"

test -r veto-definitions || git clone --depth 1 git@code.pycbc.phy.syr.edu:detchar/veto-definitions.git
VETO_PATH="file://`pwd`/veto-definitions"

echo -e "\\n>> [`date`] Building test workflow $WORKFLOWNAME"
cp `which ligo-proxy-init` .
patch -p0 ligo-proxy-init 1>/dev/null <<EOF
--- /bin/ligo-proxy-init	2016-12-05 07:18:14.000000000 -0500
+++ ligo-proxy-init	2017-04-09 12:49:35.575182509 -0400
@@ -210,7 +210,7 @@
     fi
 
     login=\${1/@*/}
-    curl_auth_method="--user \$login"
+    curl_auth_method="--user \$login:\${LIGO_TOKEN}"
     echo "Your identity: \$login@LIGO.ORG"
 fi
EOF

unset X509_USER_PROXY
export LIGO_TOKEN=`cat ~/.ssh/ldg_token`
LIGO_USER=`cat ~/.ssh/ldg_user`
./ligo-proxy-init -p ${LIGO_USER} 1>/dev/null
unset LIGO_TOKEN LIGO_USER

UUID=`uuidgen`
WORKFLOW_NAME=test-workflow-$UUID
OUTPUT_PATH=`pwd`/public_html/test_workflow/${WORKFLOW_NAME}
export LIGO_DATAFIND_SERVER="datafind.ligo.org:443"

mkdir $WORKFLOW_NAME
pushd $WORKFLOW_NAME

pycbc_make_coinc_search_workflow \
--workflow-name ${WORKFLOW_NAME} --output-dir output \
--config-files \
  ${CONFIG_PATH}/O2/pipeline/analysis.ini \
  ${CONFIG_PATH}/O2/pipeline/executables.ini \
  ${CONFIG_PATH}/O2/pipeline/injections.ini \
  ${CONFIG_PATH}/O2/pipeline/plotting.ini \
--config-overrides \
  "workflow:datafind-check-frames-exist:warn" \
  "workflow:start-time:$((1126259462 - 1800))" \
  "workflow:end-time:$((1126259462 + 1800))" \
  "results_page:output-path:${OUTPUT_PATH}" \
  "results_page:analysis-title:'PyCBC GW150914 Test Workflow'" \
  "results_page:analysis-subtitle:'$UUID'" \
  "optimal_snr:cores:24" \
  "workflow-splittable-full_data:splittable-num-banks:20" \
  "workflow-splittable-injections:splittable-num-banks:10" \
  "workflow:h1-channel-name:H1:DCS-CALIB_STRAIN_C02" \
  "workflow:l1-channel-name:L1:DCS-CALIB_STRAIN_C02" \
  "workflow-ifos:h1:" \
  "workflow-ifos:l1:" \
  "workflow-datafind:datafind-h1-frame-type:H1_HOFT_C02" \
  "workflow-datafind:datafind-l1-frame-type:L1_HOFT_C02" \
  "workflow-segments:segments-h1-science-name:H1:DCS-ANALYSIS_READY_C02:1" \
  "workflow-segments:segments-l1-science-name:L1:DCS-ANALYSIS_READY_C02:1" \
  "workflow-segments:segments-database-url:https://segments.ligo.org" \
  "workflow-segments:segments-science-veto:1" \
  "workflow-segments:segments-final-veto-group:12H" \
  "workflow-segments:segments-veto-groups:" \
  "datafind:urltype:file" \
  "workflow-segments:segments-veto-definer-url:${VETO_PATH}/cbc/O1/H1L1-CBC_VETO_DEFINER_C02_O1_1126051217-11203200.xml" \
  "workflow-tmpltbank:tmpltbank-pregenerated-bank:${CONFIG_PATH}/O1/bank/H1L1-UBERBANK_MAXM100_NS0p05_ER8HMPSD-1126033217-223200.xml.gz" \
  "inspiral:low-frequency-cutoff:30" \
  "s-mchirp:bins:0.8 1.74 8.07 14.92 21.77 100" \
  "s-mtotal:bins:2 4 27.25 51.5 75.75 100" \
--config-delete \
  "injections-bnsstt2_inj" \
  "injections-nsbhseobnrv4_inj" \
  "injections-nsbhimrphenomd_inj" \
  "injections-nsbhstt4_inj" \
  "injections-nsbhseobnrv3_inj" \
  "injections-bbhimrphenomd_inj" \
  "injections-bbhseobnrv3_inj" \
  "injections-imbhseobnrv4_inj" \
  "injections-imbhimrphenomd_inj" \
  "injections-imbhseobnrv3_inj" \
  "injections-imbheobnrv2hm_inj" \
  "inspiral:enable-bank-start-frequency"

pushd output

for workflow in *.dax
do
  echo -e "\\n>> [`date`] Validating workflow $workflow"
  pegasus-dax-validator $workflow
done

echo -e "\\n>> [`date`] Planning test workflow"
pycbc_submit_dax \
  --accounting-group ligo.dev.o2.cbc.bbh.pycbcoffline \
  --dax ${WORKFLOW_NAME}.dax \
  --no-create-proxy \
  --no-submit
popd
popd

grid-proxy-destroy

echo -e "\\n>> [`date`] Test workflow validation complete"

exit 0
