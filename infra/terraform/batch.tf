resource "aws_cloudwatch_log_group" "batch" {
  name              = "/aws/batch/${var.project_name}"
  retention_in_days = 30
}

# ── CPU compute environment ───────────────────────────────────────────────────

resource "aws_batch_compute_environment" "cpu" {
  compute_environment_name = "${var.project_name}-cpu"
  type                     = "MANAGED"
  service_role             = aws_iam_role.batch_service.arn

  compute_resources {
    type                = "SPOT"
    allocation_strategy = "SPOT_CAPACITY_OPTIMIZED"
    instance_role       = aws_iam_instance_profile.ec2_instance.arn
    instance_type       = ["m5", "c5", "r5"]
    min_vcpus           = 0
    max_vcpus           = var.max_vcpus
    subnets             = data.aws_subnets.default.ids
    security_group_ids  = [aws_security_group.batch_jobs.id]
  }
}

resource "aws_batch_job_queue" "cpu" {
  name                 = "${var.project_name}-cpu-queue"
  state                = "ENABLED"
  priority             = 1
  compute_environments = [aws_batch_compute_environment.cpu.arn]
}

resource "aws_batch_job_definition" "cpu" {
  name = "${var.project_name}-cpu-job"
  type = "container"

  container_properties = jsonencode({
    image      = "${aws_ecr_repository.compute.repository_url}:cpu"
    jobRoleArn = aws_iam_role.job.arn
    resourceRequirements = [
      { type = "VCPU",   value = "4"    },
      { type = "MEMORY", value = "8192" }
    ]
    logConfiguration = {
      logDriver = "awslogs"
      options = {
        "awslogs-group"  = aws_cloudwatch_log_group.batch.name
        "awslogs-region" = var.region
      }
    }
  })
}

# ── GPU compute environment ───────────────────────────────────────────────────

resource "aws_batch_compute_environment" "gpu" {
  compute_environment_name = "${var.project_name}-gpu"
  type                     = "MANAGED"
  service_role             = aws_iam_role.batch_service.arn

  compute_resources {
    type                = "SPOT"
    allocation_strategy = "SPOT_CAPACITY_OPTIMIZED"
    instance_role       = aws_iam_instance_profile.ec2_instance.arn
    instance_type       = ["g4dn.xlarge", "g4dn.2xlarge"]
    min_vcpus           = 0
    max_vcpus           = var.max_vcpus_gpu
    subnets             = data.aws_subnets.default.ids
    security_group_ids  = [aws_security_group.batch_jobs.id]
  }
}

resource "aws_batch_job_queue" "gpu" {
  name                 = "${var.project_name}-gpu-queue"
  state                = "ENABLED"
  priority             = 1
  compute_environments = [aws_batch_compute_environment.gpu.arn]
}

resource "aws_batch_job_definition" "gpu" {
  name = "${var.project_name}-gpu-job"
  type = "container"

  container_properties = jsonencode({
    image      = "${aws_ecr_repository.compute.repository_url}:gpu"
    jobRoleArn = aws_iam_role.job.arn
    resourceRequirements = [
      { type = "VCPU",   value = "8"     },
      { type = "MEMORY", value = "32768" },
      { type = "GPU",    value = "1"     }
    ]
    logConfiguration = {
      logDriver = "awslogs"
      options = {
        "awslogs-group"  = aws_cloudwatch_log_group.batch.name
        "awslogs-region" = var.region
      }
    }
  })
}
