resource "aws_ecr_repository" "compute" {
  name                 = "${var.project_name}-compute"
  image_tag_mutability = "MUTABLE"

  image_scanning_configuration {
    scan_on_push = true
  }
}

resource "aws_ecr_lifecycle_policy" "compute" {
  repository = aws_ecr_repository.compute.name

  policy = jsonencode({
    rules = [
      {
        rulePriority = 1
        description  = "Keep last 5 images per tag prefix"
        selection = {
          tagStatus     = "tagged"
          tagPrefixList = ["cpu", "gpu"]
          countType     = "imageCountMoreThan"
          countNumber   = 5
        }
        action = { type = "expire" }
      }
    ]
  })
}
