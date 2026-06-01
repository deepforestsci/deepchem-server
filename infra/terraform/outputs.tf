output "bucket_name" {
  description = "S3 data bucket name — set as AWS_BUCKET"
  value       = aws_s3_bucket.data.bucket
}

output "ecr_url" {
  description = "ECR repository URL — used by build_and_push.sh"
  value       = aws_ecr_repository.compute.repository_url
}

output "cpu_job_queue" {
  description = "CPU Batch job queue name — set as AWS_BATCH_CPU_JOB_QUEUE"
  value       = aws_batch_job_queue.cpu.name
}

output "gpu_job_queue" {
  description = "GPU Batch job queue name — set as AWS_BATCH_GPU_JOB_QUEUE"
  value       = aws_batch_job_queue.gpu.name
}

output "cpu_job_definition" {
  description = "CPU Batch job definition name — set as AWS_BATCH_CPU_JOB_DEFINITION"
  value       = aws_batch_job_definition.cpu.name
}

output "gpu_job_definition" {
  description = "GPU Batch job definition name — set as AWS_BATCH_GPU_JOB_DEFINITION"
  value       = aws_batch_job_definition.gpu.name
}
