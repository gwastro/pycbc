#!/bin/bash

set -e

echo -e "\\n>> [`date`] Testing pycbc_make_coinc_search_workflow"

VENV_PATH=${1}
TRAVIS_TAG=${2}

if [ "x${VENV_PATH}" == "x" ] ; then
  echo -e "\\n>> [`date`] Error: VENV_PATH was not passed to script or is empty"
  exit 1
fi

if [ "x${TRAVIS_TAG}" == "x" ] ; then
  echo -e "\\n>> [`date`] Error: TRAVIS_TAG was not passed to script or is empty"
  exit 1
fi

echo -e "\\n>> [`date`] Entering virtual environment $VENV_PATH"
source ${VENV_PATH}/bin/activate

CONFIG_PATH="https://raw.githubusercontent.com/ligo-cbc/pycbc-config/${TRAVIS_TAG}"
echo -e "\\n>> [`date`] Using config files from ${CONFIG_PATH}"

VETO_DEFINER="https://git.ligo.org/detchar/veto-definitions/raw/a07f542b37ccfcfbf2b732f0d75d0f1ab4166d9f/cbc/O1/H1L1-CBC_VETO_DEFINER_C02_O1_1126051217-11203200.xml"
echo -e "\\n>> [`date`] Using veto definer file from ${VETO_DEFINER}"

BANK_FILE="${CONFIG_PATH}/O1/bank/H1L1-UBERBANK_MAXM100_NS0p05_ER8HMPSD-1126033217-223200.xml.gz"
echo -e "\\n>> [`date`] Using template bank from ${BANK_FILE}"

echo -e "\\n>> [`date`] Patching ligo-proxy-init for Travis"
cp `which ligo-proxy-init` .
set +e
patch -p0 ligo-proxy-init <<EOF
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
if [ ! $? -eq 0 ] ; then
patch -p0 ligo-proxy-init <<EOF
--- /bin/ligo-proxy-init	2017-04-12 12:27:45.000000000 +0000
+++ ligo-proxy-init	2017-09-07 23:37:51.224116188 +0000
@@ -212,7 +212,7 @@
 
     login=\${1/@*/}
     [[ \$login == *","* ]] && echo "Replacing comma characters in login!"; login=\${login//,/.}
-    curl_auth_method="--user \$login"
+    curl_auth_method="--user \$login:\${LIGO_TOKEN}"
     echo "Your identity: \$login@LIGO.ORG"
 fi
EOF
if [ ! $? -eq 0 ] ; then
echo -e "\\n>> [`date`] ERROR: could not patch ligo-proxy-init for Travis"
exit 1
fi
fi
set -e

echo -e "\\n>> [`date`] Patching ecp-cookie-init for Travis"
cp `which ecp-cookie-init` .
patch -p0 ecp-cookie-init << EOF
--- /bin/ecp-cookie-init        2016-12-21 08:41:13.000000000 -0500
+++ /tmp/ecp-cookie-init        2017-07-11 09:43:30.846451317 -0400
@@ -268,7 +268,7 @@
     target=\$2
     login=\$3
 
-    curl_auth_method="--user \$login"
+    curl_auth_method="--user \$login:\${LIGO_TOKEN}"
 fi
 
 if [ "\${idp_tag}" = "LIGO.ORG" ] ; then
EOF

unset X509_USER_PROXY
export LIGO_TOKEN=`cat ~/.ssh/ldg_token`
LIGO_USER=`cat ~/.ssh/ldg_user`

echo -e "\\n>> [`date`] Creating proxy"
./ligo-proxy-init -p ${LIGO_USER} 1>/dev/null

echo -e "\\n>> [`date`] Creating ECP cookie"
./ecp-cookie-init LIGO.ORG https://git.ligo.org/users/auth/shibboleth/callback ${LIGO_USER} 1>/dev/null
unset LIGO_TOKEN LIGO_USER

echo -e "\\n>> [`date`] Creating test workflow"
UUID=`uuidgen`
WORKFLOW_NAME=test-workflow-$UUID
OUTPUT_PATH=`pwd`/public_html/test_workflow/${WORKFLOW_NAME}
export LIGO_DATAFIND_SERVER="datafind.ligo.org:443"

mkdir $WORKFLOW_NAME
pushd $WORKFLOW_NAME

echo -e "\\n>> [`date`] Building test workflow $WORKFLOWNAME"

set +e

pycbc_make_coinc_search_workflow \
--workflow-name ${WORKFLOW_NAME} --output-dir output \
--config-files \
  ${CONFIG_PATH}/O2/pipeline/analysis.ini \
  ${CONFIG_PATH}/O2/pipeline/executables.ini \
  ${CONFIG_PATH}/O2/pipeline/injections.ini \
  ${CONFIG_PATH}/O2/pipeline/plotting.ini \
--config-overrides \
  "workflow-datafind:datafind-check-frames-exist:warn" \
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
  "workflow-segments:segments-veto-definer-url:${VETO_DEFINER}" \
  "workflow-tmpltbank:tmpltbank-pregenerated-bank:${BANK_FILE}" \
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

if [ -f /pycbc/workflow-test/test-workflow-*/output/results/1._analysis_time/1.01_segment_data/logs/segments_from_cats.err ] ; then
  cat /pycbc/workflow-test/test-workflow-*/output/results/1._analysis_time/1.01_segment_data/logs/segments_from_cats.err
  cat /pycbc/workflow-test/test-workflow-*/output/results/1._analysis_time/1.01_segment_data/logs/segments_from_cats.out
fi

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

set -e 

grid-proxy-destroy

echo -e "\\n>> [`date`] Test workflow validation complete"

exit 0
