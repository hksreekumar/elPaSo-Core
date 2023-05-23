# Container Registry

Container registries on the elPaSo GitLab repository hold basic OS and CI deployed elPaSo docker images.
 - **OS images** supply the desired environment to build and run elPaSo, whereas
 - **elPaSo images** contain CI built images with elPaSo ready to use.

(docker-setup)=
## Docker setup

{bdg-secondary}`stable | Docker 20.10.21`

Refer to [Docker Documentation](https://docs.docker.com/) for detailed understanding. Follow the steps below for the initial setup:
- See installation requirements and instructions from docker [here](https://docs.docker.com/engine/install/ubuntu/).
    ```bash
    # remove older versions
    sudo apt-get remove docker docker-engine docker.io containerd runc
    #
    sudo apt-get upgrade
    sudo apt-get update
    sudo apt-get install ca-certificates curl gnupg lsb-release
    sudo mkdir -p /etc/apt/keyrings
    sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /etc/apt/keyrings/docker.gpg
    sudo echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu  $(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null
    sudo apt-get update
    sudo chmod a+r /etc/apt/keyrings/docker.gpg
    # 
    sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
    sudo docker run hello-world
    ```
- If you desire to push new containers into the registry, you require access token generated for role `developer` with scope `api`:
    ```bash
    sudo docker login -u gitlab-ci-token -p $TOKEN git.rz.tu-bs.de:4567
    ```

### Docker cheatsheet

Following is a consolidated list of useful Docker commands:
```bash
# list all available local images
docker images

# list all available local containers
docker ps

# build an image
docker build . -t name_of_docker_image:version

# pull and push an image from and to a registry (Default is to look in DockerHub)
docker pull name_of_docker_image:version
docker push name_of_docker_image:version

# hosting container in GitLab's container registry
docker login git.rz.tu-bs.de:4567 -p <ACCESS_TOKEN>
docker build . -t git.rz.tu-bs.de:4567/<PROJECT_PATH>:VERSION
docker push git.rz.tu-bs.de:4567/<PROJECT_PATH>:VERSION

# run a container using interactive bash shell
docker run -it name_of_docker_image:version

# copy a file from local system to a container
docker cp <SRC_PATH_TO_FILE_IN_SYSTEM> <CONTAINER_ID>:<DESTINATION_PATH_IN_CONTAINER>

# stop a container
docker stop <CONTAINER_ID>

# save and load docker image as tar
docker save --output archive.tar <IMAGE>
docker load --input archive.tar

# cleaning
sudo docker image prune -a
sudo docker system prune
sudo docker system df
```

## Maintaining elPaSo containers

### OS images

```bash
cd <ELPASO_REPOSITORY_PATH>/.devcontainer/docker-gnu
sudo docker build -t git.rz.tu-bs.de:4567/akustik/elpaso-core/elpasocore-baseimage-ubuntu-x64:23.01.1 .
sudo docker push git.rz.tu-bs.de:4567/akustik/elpaso-core/elpasocore-baseimage-ubuntu-x64:23.01.1
```

### elPaSo images

Elpaso images are built by are deployed by the CI pipeline after sucessful run. Image name is: `git.rz.tu-bs.de:4567/akustik/elpaso-core/elpasocore:latest`.

Use the latest image with:
```bash
sudo docker run -it git.rz.tu-bs.de:4567/akustik/elpaso-core/elpasocore:latest
```