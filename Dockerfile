FROM quay.io/biocontainers/hapcut2:1.3.3--hb0d9459_3
WORKDIR /usr/local/hapcut2_mec_solver
RUN pip install pip==22.3.1 setuptools==67.0.0 pandas==1.5.3 pytest==7.2.1 scipy==1.10.0
COPY . /usr/local/hapcut2_mec_solver
RUN pip install -v .
RUN pytest
RUN python -m hapcut2_mec_solver < tests/test.json


