FROM pycbc/pycbc-base-el7:v1.9-89a30fc

RUN echo '#!/bin/bash -l' > /opt/pycbc/entrypoint.sh
RUN echo 'exec $@' >> /opt/pycbc/entrypoint.sh
RUN chmod 0755 /opt/pycbc/entrypoint.sh
ENTRYPOINT ["/opt/pycbc/entrypoint.sh"]

USER pycbc
WORKDIR /opt/pycbc
