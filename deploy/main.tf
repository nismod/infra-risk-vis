//
// Terraform script - only used in AWS deployment.
//


//
// Variables
//

variable "SITE_URL" {
  default = "jamaica.infrastructureresilience.org"
}

variable "instance_type" {
  description = "AWS EC2 instance type"
  default = "t3.micro"
}

provider "aws" {
  region = "eu-west-2"
}

//
// EC2 Connection
// Keypair, VPC, Security Group
//

resource "aws_key_pair" "deployer" {
  key_name   = "opsis-aws-deployer-jamaica"
  public_key = "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIAjX5yhc7GWROTVOM8r92rO6MEUyKt/JfTCQzCY/lNi9 opsis-aws-jamaica"
}

resource "aws_default_vpc" "jamaica" {
  tags = {
    Name = "Jamaica VPC"
  }
}

resource "aws_security_group" "access_http_ssh" {
  name = "access_jamaica"
  vpc_id = aws_default_vpc.jamaica.id
  ingress {
    from_port   = 22
    to_port     = 22
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }
  ingress {
    from_port   = 443
    to_port     = 443
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }
  ingress {
    from_port   = 80
    to_port     = 80
    protocol    = "tcp"

    cidr_blocks = ["0.0.0.0/0"]
  }
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

//
// EC2 Instance
//

data "aws_ami" "ubuntu" {
  most_recent = true

  filter {
    name   = "name"
    values = ["ubuntu/images/hvm-ssd/ubuntu-focal-20.04-amd64-server-*"]
  }

  filter {
    name   = "virtualization-type"
    values = ["hvm"]
  }

  owners = ["099720109477"] # Canonical
}

resource "aws_instance" "jamaica" {
  instance_type = var.instance_type
  ami = data.aws_ami.ubuntu.id

  key_name = aws_key_pair.deployer.key_name
  vpc_security_group_ids = [aws_security_group.access_http_ssh.id]

  tags = {
    Name = "${var.SITE_URL} ${var.instance_type}"
  }
}

//
// DNS
//

data "aws_route53_zone" "selected" {
  name = "infrastructureresilience.org."
}

resource "aws_route53_record" "jamaica" {
  zone_id = data.aws_route53_zone.selected.zone_id
  name    = var.SITE_URL
  type    = "A"
  ttl     = "300"
  records = [aws_instance.jamaica.public_ip]
}

output "jamaica-public_ip" {
  value = aws_instance.jamaica.public_ip
}
