# SjS-DigitalTwin — Image reproductible
#
# Base : colomoto/colomoto-docker (inclut CaSQ, MaBoSS, GINsim, bioLQM, pyboolnet, mpbn)
# Couche projet : env conda sjs-digitaltwin (envs/environment.yml)
#
# Build :
#   docker build -t sjs-digitaltwin:latest .
#
# Run interactif :
#   docker run --rm -it -v "$PWD":/work -w /work sjs-digitaltwin:latest bash
#
# Pipeline non-interactif :
#   docker run --rm -v "$PWD":/work -w /work sjs-digitaltwin:latest \
#       conda run -n sjs-digitaltwin python scripts/01_download_sjd_map.py

FROM colomoto/colomoto-docker:2024-08-01

LABEL org.opencontainers.image.title="SjS-DigitalTwin"
LABEL org.opencontainers.image.description="Boolean immune digital twin of the salivary gland in Sjogren's disease"
LABEL org.opencontainers.image.source="https://github.com/Nurtal/IDT_SjS_gland"
LABEL org.opencontainers.image.licenses="MIT"

USER root

# Outils système utiles (curl pour HTTP debug, git pour versioning, build-essential pour quelques pip wheels)
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        git \
        git-lfs \
        build-essential \
    && git lfs install \
    && rm -rf /var/lib/apt/lists/*

# Création de l'environnement conda projet
COPY envs/environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && \
    mamba clean -afy && \
    rm /tmp/environment.yml

# Activation par défaut de l'env projet pour les sessions interactives
RUN echo "conda activate sjs-digitaltwin" >> /etc/skel/.bashrc && \
    echo "conda activate sjs-digitaltwin" >> /root/.bashrc

# Variables d'environnement
ENV CONDA_DEFAULT_ENV=sjs-digitaltwin
ENV PATH=/opt/conda/envs/sjs-digitaltwin/bin:$PATH
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

WORKDIR /work

# Healthcheck minimal : import des modules clés
HEALTHCHECK --interval=60s --timeout=10s --start-period=5s --retries=3 \
    CMD conda run -n sjs-digitaltwin python -c "import casq, libsbml, pyboolnet, mpbn" || exit 1

CMD ["bash"]
