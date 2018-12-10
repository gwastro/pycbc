FROM pycbc/pycbc-base-el7:v3.2-89a30fc

ENTRYPOINT ["/bin/bash", "-vc", "source /opt/pycbc/pycbc-software/bin/activate; exec $0 $@"]

RUN echo '#!/bin/bash -l' > /opt/pycbc/entrypoint.sh
RUN echo 'exec $@' >> /opt/pycbc/entrypoint.sh
RUN chmod 0755 /opt/pycbc/entrypoint.sh
ENTRYPOINT ["/opt/pycbc/entrypoint.sh"]

USER pycbc
WORKDIR /opt/pycbc
