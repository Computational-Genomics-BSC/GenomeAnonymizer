# Multi-stage build for GenomeAnonymizer
# Stage 1: Build the application
FROM gradle:8.5-jdk21 AS builder

# Set working directory
WORKDIR /app

# Copy gradle files first for better layer caching
COPY build.gradle.kts ./
COPY gradle/ ./gradle/
COPY gradlew ./
COPY gradlew.bat ./

# Make gradlew executable
RUN chmod +x ./gradlew

# Copy source code
COPY src/ ./src/

# Build the application using gradlew (creates fat JAR with all dependencies)
RUN ./gradlew clean build --no-daemon --info

# List the contents of build/libs to debug
RUN ls -la /app/build/libs/

# Stage 2: Runtime image
FROM eclipse-temurin:21-jre

# Set working directory
WORKDIR /app

# Copy the built JAR from builder stage
COPY --from=builder /app/build/libs/app*.jar /GenomeAnonymizer.jar

# Create an entrypoint script
#RUN echo '#!/bin/bash\njava "$@" -jar /app/GenomeAnonymizer.jar' > /app/entrypoint.sh && \
#    chmod +x /app/entrypoint.sh

# Set the entrypoint
#ENTRYPOINT ["/app/entrypoint.sh"]

# Default command shows help
#CMD ["--help"]
CMD [ "sh" ]

# Add labels for metadata
LABEL maintainer="Nicolas Gaitan & Rodrigo Martin - BSC-CNS Computational Genomics Group"
LABEL description="GenomeAnonymizer - Software for anonymizing Short Read WGS/WES data"
LABEL version="1.1.0"
LABEL java.version="21"
