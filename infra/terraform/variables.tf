variable "region" {
  description = "AWS region"
  default     = "us-east-1"
}

variable "project_name" {
  description = "Prefix for all resource names"
  default     = "deepchem-server"
}

variable "bucket_name" {
  description = "Override S3 bucket name. Leave empty to use account-ID-based default."
  default     = ""
}

variable "max_vcpus" {
  description = "Max vCPUs for CPU compute environment"
  default     = 256
}

variable "max_vcpus_gpu" {
  description = "Max vCPUs for GPU compute environment"
  default     = 64
}
