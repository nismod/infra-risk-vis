//
// Versions for terraform plugins
//
terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = ">= 3.20.0"
    }
  }
}


//
// Variables
//
variable "SITE_URL" {
  default = "caribbean.infrastructureresilience.org"
}

variable "instance_type" {
  description = "AWS EC2 instance type"
  default     = "t3.micro"
}

provider "aws" {
  region = "eu-west-2"
}

// Leave this empty, to be prompted when running terraform apply
variable "RDS_PASSWORD" {}

//
// EC2 Connection
// Keypair, VPC, Security Group
//

resource "aws_key_pair" "deployer" {
  key_name   = "opsis-aws-deployer-caribbean"
  public_key = "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIBhi88MtbHeh0FwIQA4bb6DOS5SQZ61buTooGPOUfXv7 opsis-aws-caribbean"
}

data "aws_availability_zones" "available" {}

module "vpc" {
  source  = "terraform-aws-modules/vpc/aws"
  version = "2.77.0"

  name                 = "caribbean"
  cidr                 = "10.0.0.0/16"
  azs                  = data.aws_availability_zones.available.names
  public_subnets       = ["10.0.4.0/24", "10.0.5.0/24", "10.0.6.0/24"]
  enable_dns_hostnames = true
  enable_dns_support   = true
}

resource "aws_db_subnet_group" "caribbean" {
  name       = "caribbean"
  subnet_ids = module.vpc.public_subnets

  tags = {
    Name = "caribbean"
  }
}

// CIDR blocks below are too permissive - edit to constrain to specific IPs
// or ranges.
// - web access from anywhere
// - SSH access from dev machines if at all
// - DB access from VM (within subnet) and possibly dev machines
resource "aws_security_group" "access" {
  name   = "access_caribbean"
  vpc_id = module.vpc.vpc_id
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
  ingress {
    from_port   = 5432
    to_port     = 5432
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }

  egress {
    from_port   = 5432
    to_port     = 5432
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

resource "aws_instance" "caribbean" {
  instance_type = var.instance_type
  ami = data.aws_ami.ubuntu.id

  key_name = aws_key_pair.deployer.key_name
  subnet_id = module.vpc.public_subnets[0]
  vpc_security_group_ids = [aws_security_group.access.id]

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

resource "aws_route53_record" "caribbean" {
  zone_id = data.aws_route53_zone.selected.zone_id
  name    = var.SITE_URL
  type    = "A"
  ttl     = "300"
  records = [aws_instance.caribbean.public_ip]
}

output "caribbean-public_ip" {
  value = aws_instance.caribbean.public_ip
}

//
// Database
//

resource "aws_db_parameter_group" "pg-caribbean-dev" {
  name   = "pg-caribbean-dev"
  family = "postgres14"

  parameter {
    name  = "log_connections"
    value = "1"
  }
}

resource "aws_db_instance" "pg-caribbean-dev" {
  allocated_storage     = 20
  max_allocated_storage = 100

  engine               = "postgres"
  engine_version       = "14.2"

  instance_class = "db.t3.micro"
  identifier     = "pg-caribbean-dev"
  name           = "caribbean_dev"
  username       = "caribbean_dev"
  port           = 5432
  password       = var.RDS_PASSWORD

  multi_az               = false

  db_subnet_group_name   = aws_db_subnet_group.caribbean.name
  vpc_security_group_ids = [aws_security_group.access.id]
  parameter_group_name   = aws_db_parameter_group.pg-caribbean-dev.name
  publicly_accessible    = true
  skip_final_snapshot    = true

  backup_retention_period = 2
  tags = {
    Name = "caribbean_dev-pg-instance"
  }
}

output "pg-caribbean-dev_hostname" {
  description = "RDS instance hostname"
  value       = aws_db_instance.pg-caribbean-dev.address
  sensitive   = true
}

output "pg-caribbean-dev_port" {
  description = "RDS instance port"
  value       = aws_db_instance.pg-caribbean-dev.port
  sensitive   = true
}

output "pg-caribbean-dev_username" {
  description = "RDS instance root username"
  value       = aws_db_instance.pg-caribbean-dev.username
  sensitive   = true
}
