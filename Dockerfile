FROM quay.io/biocontainers/hapcut2:1.3.3--hb0d9459_3
COPY . /usr/local/hapcut2_mec_solver
WORKDIR /usr/local/hapcut2_mec_solver
RUN pip install pip==22.3.1 setuptools==67.0.0
RUN pip install -v .
RUN chmod +x tests/test.sh


