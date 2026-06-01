output "bucket_name" {
  description = "S3 data bucket name — set as AWS_BUCKET env var"
  value       = aws_s3_bucket.data.bucket
}

output "ecr_url" {
  description = "ECR repository URL — used by build_and_push.sh"
  value       = aws_ecr_repository.compute.repository_url
}
