FROM pycbc/pycbc-base-el7:v1.9-89a30fc

USER pycbc
WORKDIR /opt/pycbc

RUN echo '#!/bin/bash -l' > /etc/entrypoint.sh
RUN echo 'exec $@' >> /etc/entrypoint.sh
RUN chmod 0755 /etc/entrypoint.sh
ENTRYPOINT ["/etc/entrypoint.sh"]
