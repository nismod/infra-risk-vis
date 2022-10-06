//
// Variables
//

variable "SITE_URL" {
  default = "global.infrastructureresilience.org"
}

variable "instance_type" {
  description = "AWS EC2 instance type"
  default     = "t3.micro"
}

provider "aws" {
  region = "eu-west-2"
}

// Leave this empty, to be prompted when running terraform apply
variable "RDS_POSTGRES_PASSWORD" {}
variable "RDS_MYSQL_PASSWORD" {}

//
// EC2 Connection
// Keypair, VPC, Security Group
//

resource "aws_key_pair" "deployer" {
  key_name   = "opsis-aws-deployer-global"
  public_key = "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIJi8SKyeEXTv1Mev8BU4iH6xGb3PRyNnDtZbbjQcsHLp opsis-aws-deployer-global"
}

data "aws_availability_zones" "available" {}

module "vpc" {
  source  = "terraform-aws-modules/vpc/aws"
  version = "2.77.0"

  name                 = "global"
  cidr                 = "10.0.0.0/16"
  azs                  = data.aws_availability_zones.available.names
  public_subnets       = ["10.0.4.0/24", "10.0.5.0/24", "10.0.6.0/24"]
  enable_dns_hostnames = true
  enable_dns_support   = true
}

resource "aws_db_subnet_group" "global" {
  name       = "global"
  subnet_ids = module.vpc.public_subnets

  tags = {
    Name = "Global VPC"
  }
}

resource "aws_security_group" "access_http_ssh" {
  name   = "access_global"
  vpc_id = "vpc-b39bfcdb"
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
    from_port = 80
    to_port   = 80
    protocol  = "tcp"

    cidr_blocks = ["0.0.0.0/0"]
  }
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

resource "aws_security_group" "access" {
  name   = "access_global"
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
    from_port = 80
    to_port   = 80
    protocol  = "tcp"

    cidr_blocks = ["0.0.0.0/0"]
  }
  ingress {
    from_port   = 5432
    to_port     = 5432
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }
  ingress {
    from_port   = 3306
    to_port     = 3306
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
  egress {
    from_port   = 5432
    to_port     = 5432
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }
  egress {
    from_port   = 3306
    to_port     = 3306
    protocol    = "tcp"
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

resource "aws_instance" "global" {
  instance_type = var.instance_type
  ami           = data.aws_ami.ubuntu.id

  key_name               = aws_key_pair.deployer.key_name
  subnet_id              = module.vpc.public_subnets[0]
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

resource "aws_route53_record" "global" {
  zone_id = data.aws_route53_zone.selected.zone_id
  name    = var.SITE_URL
  type    = "A"
  ttl     = "300"
  records = [aws_instance.global.public_ip]
}

output "global-public_ip" {
  value = aws_instance.global.public_ip
}

//
// Postgres Database
//

resource "aws_db_parameter_group" "pg-global-dev" {
  name   = "pg-global-dev"
  family = "postgres14"

  parameter {
    name  = "log_connections"
    value = "1"
  }
}

resource "aws_db_instance" "pg-global-dev" {
  allocated_storage     = 20
  max_allocated_storage = 200

  engine         = "postgres"
  engine_version = "14.2"

  instance_class = "db.t3.micro"
  identifier     = "pg-global-dev"
  db_name        = "global_dev"
  username       = "global_dev"
  port           = 5432
  password       = var.RDS_POSTGRES_PASSWORD

  multi_az = false

  db_subnet_group_name   = aws_db_subnet_group.global.name
  vpc_security_group_ids = [aws_security_group.access.id]
  parameter_group_name   = aws_db_parameter_group.pg-global-dev.name
  publicly_accessible    = true
  skip_final_snapshot    = true

  backup_retention_period = 2
  tags = {
    Name = "global_dev-pg-instance"
  }
}

output "pg-global-dev_hostname" {
  description = "PG RDS instance hostname"
  value       = aws_db_instance.pg-global-dev.address
  sensitive   = true
}

output "pg-global-dev_port" {
  description = "PG RDS instance port"
  value       = aws_db_instance.pg-global-dev.port
  sensitive   = true
}

output "pg-global-dev_username" {
  description = "PG RDS instance root username"
  value       = aws_db_instance.pg-global-dev.username
  sensitive   = true
}

//
// MySQL Database
//

resource "aws_db_parameter_group" "mysql-global-dev" {
  name   = "mysql-global-dev"
  family = "mysql8.0"

}

resource "aws_db_instance" "mysql-global-dev" {
  allocated_storage     = 20
  max_allocated_storage = 100

  engine         = "mysql"
  engine_version = "8.0.27"

  instance_class = "db.t3.micro"
  identifier     = "mysql-global-dev"
  db_name        = "global_dev"
  username       = "global_dev"
  port           = 3306
  password       = var.RDS_MYSQL_PASSWORD

  multi_az = false

  db_subnet_group_name   = aws_db_subnet_group.global.name
  vpc_security_group_ids = [aws_security_group.access.id]
  parameter_group_name   = aws_db_parameter_group.mysql-global-dev.name
  publicly_accessible    = true
  skip_final_snapshot    = true

  backup_retention_period = 2
  tags = {
    Name = "global_dev-mysql-instance"
  }
}

output "mysql-global-dev_hostname" {
  description = "MYSQL RDS instance hostname"
  value       = aws_db_instance.mysql-global-dev.address
  sensitive   = true
}

output "mysql-global-dev_port" {
  description = "MYSQL RDS instance port"
  value       = aws_db_instance.mysql-global-dev.port
  sensitive   = true
}

output "mysql-global-dev_username" {
  description = "MYSQL RDS instance root username"
  value       = aws_db_instance.mysql-global-dev.username
  sensitive   = true
}
