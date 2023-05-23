# CI Pipeline

Refer to [GitLab CI Documentation](https://docs.gitlab.com/ee/ci/) for basic understanding.

elPaSo has a classical CI pipeline that performs automatic building, testing and deployment. See pipeline in **elPaSo repository > CI/CD > Pipelines**. Gitlab `yaml` scripts also serve as exact  readable documentation of the entire pipeline. Refer to elPaSo CI-yaml scripts: `.gitlab-ci.yml, .gitlab-ci-pipeline-intel.yml, .gitlab-ci-pipeline-gnu.yml`.

## GitLab runner configuration

{bdg-secondary}`stable | gitlab-runner 14.8.2`

Further details in [Gitlab runner documentation](https://docs.gitlab.com/runner/).

elPaSo currently perform its CI/CD pipeline through the following runner:
- Institute workstation; #DockerExecutor #Ubuntu. Following are the steps to set up the runner:
    - Runner installation; provide runner tags `DockerExecutor, Ubuntu`
        ```bash
        curl -L "https://packages.gitlab.com/install/repositories/runner/gitlab-runner/script.deb.sh" | sudo bash
        sudo apt-get install gitlab-runner
        sudo gitlab-runner register # follow basic instructions and add tags "DockerExecutor, Ubuntu"
        ```
    - Some additions to the runner configuration. Add the following configs to the runner config file `/etc/gitlab-runner/config.toml`
        ```config
        [[runners]]
            name = "XYZ"
            output_limit = 1000000
            url = "https://git.rz.tu-bs.de/"
            token = "XYZ"
            executor = "docker"
            [runners.custom_build_dir]
                enabled = true
            [runners.cache]
                [runners.cache.s3]
                [runners.cache.gcs]
                [runners.cache.azure]
            [runners.docker]
                tls_verify = false
                image = "git.rz.tu-bs.de:4567/akustik/elpaso-core/elpasocore-baseimage-ubuntu-x64"
                privileged = false
                disable_entrypoint_overwrite = false
                oom_kill_disable = false
                disable_cache = false
                volumes = ["/var/run/docker.sock:/var/run/docker.sock", "/cache"]
                pull_policy = ["if-not-present"]
                shm_size = 0
        ```
    - Start the runner
        ```bash
        sudo gitlab-runner start
        ```
    - Configuration in repository page. Follow **Project repository > Settings > CI/CD > Runners > Edit (the added runner) > Check "Run untagged jobs"**


## Manual job execution

If you would like to run the jobs individually on your local machine, execute:
```bash
sudo gitlab-runner exec docker my-job-name 
```

## Required variables

The following variable are required to be set that are called during the CI pipeline. In **Project repository > Settings > CI/CD > Variables**
1. <tt>CI_CONAN_MANAGER</tt> with token for conan deployment
2. <tt>CI_JOB_TOKEN</tt> with token for docker deployment
3. <tt>CI_JOB_TOKEN_AUTOMATE_ISSUE</tt> with token for bug reporting

## elPaSo image with intel compilers for CI run

A non-distributed image is loaded locally in the running using the following command:

```bash
sudo docker pull git.rz.tu-bs.de:4567/akustik/elpaso-research/elpaso-intelstudio-ubuntu-x64:23.01.1
sudo docker tag git.rz.tu-bs.de:4567/akustik/elpaso-research/elpaso-intelstudio-ubuntu-x64:23.01.1 localhost:5000/elpaso-intel-ubuntu-x64:22.12.1
```