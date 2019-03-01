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

VETO_DEFINER="https://raw.githubusercontent.com/gwastro/pycbc-config/master/O1/dq/H1L1-DUMMY_O1_CBC_VDEF-1126051217-1220400.xml"
echo -e "\\n>> [`date`] Using veto definer file from ${VETO_DEFINER}"

BANK_FILE="${CONFIG_PATH}/O1/bank/H1L1-UBERBANK_MAXM100_NS0p05_ER8HMPSD-1126033217-223200.xml.gz"
echo -e "\\n>> [`date`] Using template bank from ${BANK_FILE}"

LOSC_SEG_QUERY=`which pycbc_losc_segment_query`

echo -e "\\n>> [`date`] Creating test workflow"
UUID=`uuidgen`
WORKFLOW_NAME=test-workflow-$UUID
OUTPUT_PATH=`pwd`/public_html/test_workflow/${WORKFLOW_NAME}
export LIGO_DATAFIND_SERVER='128.230.190.43:80'

mkdir $WORKFLOW_NAME
pushd $WORKFLOW_NAME

echo -e "\\n>> [`date`] Building test workflow $WORKFLOWNAME"

pycbc_make_coinc_search_workflow \
--workflow-name ${WORKFLOW_NAME} --output-dir output \
--config-files \
  ${CONFIG_PATH}/O2/pipeline/analysis.ini \
  ${CONFIG_PATH}/O2/pipeline/executables.ini \
  ${CONFIG_PATH}/O2/pipeline/injections.ini \
  ${CONFIG_PATH}/O2/pipeline/plotting.ini \
--config-overrides \
  "executables:segment_query:${LOSC_SEG_QUERY}" \
  "workflow-datafind:datafind-check-frames-exist:warn" \
  "workflow:start-time:$((1126259462 - 1800))" \
  "workflow:end-time:$((1126259462 + 1800))" \
  "results_page:output-path:${OUTPUT_PATH}" \
  "results_page:analysis-title:'PyCBC GW150914 Test Workflow'" \
  "results_page:analysis-subtitle:'$UUID'" \
  "optimal_snr:cores:24" \
  "workflow-splittable-full_data:splittable-num-banks:20" \
  "workflow-splittable-injections:splittable-num-banks:10" \
  "workflow:h1-channel-name:H1:GWOSC-16KHZ_R1_STRAIN" \
  "workflow:l1-channel-name:L1:GWOSC-16KHZ_R1_STRAIN" \
  "workflow-ifos:h1:" \
  "workflow-ifos:l1:" \
  "workflow-datafind:datafind-h1-frame-type:H1_LOSC_16_V1" \
  "workflow-datafind:datafind-l1-frame-type:L1_LOSC_16_V1" \
  "workflow-segments:segments-h1-science-name:H1:RESULT:1" \
  "workflow-segments:segments-l1-science-name:L1:RESULT:1" \
  "workflow-segments:segments-database-url:https://losc.ligo.org/archive/O1/" \
  "workflow-segments:segments-science-veto:1" \
  "workflow-segments:segments-final-veto-group:12H" \
  "workflow-segments:segments-veto-groups:" \
  "workflow-segments:segments-generate-segment-files:if_not_present" \
  "datafind:urltype:file" \
  "workflow-segments:segments-veto-definer-url:${VETO_DEFINER}" \
  "workflow-tmpltbank:tmpltbank-pregenerated-bank:${BANK_FILE}" \
  "inspiral:low-frequency-cutoff:30" \
  "s-mchirp:bins:0.8 1.74 8.07 14.92 21.77 100" \
  "s-mtotal:bins:2 4 27.25 51.5 75.75 100" \
  "coinc:statistic-files:https://raw.githubusercontent.com/gwastro/pycbc-config/master/O2/pipeline/H1L1-PHASE_TIME_AMP_dummy.hdf" \
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
  --force-no-accounting-group \
  --dax ${WORKFLOW_NAME}.dax \
  --no-create-proxy \
  --no-submit \
  --no-grid

popd
popd

echo -e "\\n>> [`date`] Test workflow validation complete"

exit 0
