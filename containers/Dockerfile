FROM condaforge/miniforge3:24.3.0-0
WORKDIR /opt/plantvirome
COPY environment.yml pyproject.toml README.md LICENSE ./
COPY src ./src
RUN mamba env create -n plantvirome -f environment.yml && mamba clean -afy
ENV PATH=/opt/conda/envs/plantvirome/bin:$PATH
ENTRYPOINT ["plantvirome"]
