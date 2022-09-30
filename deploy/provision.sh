#!/usr/bin/env bash

#
# Provision virtual machine
# - assuming OS is Ubuntu 20.04 LTS
#

# Install apt packages (Docker, Compose)
sudo apt-get update
sudo apt-get install -y \
    certbot \
    ca-certificates \
    curl \
    gnupg \
    lsb-release
sudo mkdir -p /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin

# Set up SSL
# TODO - Need Oxford Uni email for reg 
sudo certbot certonly --standalone --preferred-challenges http --register-unsafely-without-email -d global.infrastructureresilience.org
# Certs arrive here: /etc/letsencrypt/live/global.infrastructureresilience.org

## The following are first-setup only - they are managed by git actions following provision
# git action user and docker group
sudo groupadd docker
sudo useradd -G docker gitaction
sudo usermod -s /bin/bash gitaction
# Generate SSH key for use with appleboy-gitaction
sudo su gitaction -c "ssh-keygen -t rsa -N '' -f ~/.ssh"
# Hosting Area
sudo mkdir /opt/infra-risk
sudo chown gitaction:gitaction /opt/infra-risk

