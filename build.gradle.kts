plugins {
    id("java")
    application
}

group = "org.example"
version = "1.1.0"

java {
    toolchain {
        languageVersion.set(JavaLanguageVersion.of(21))
    }
}

repositories {
    mavenCentral()
}

dependencies {
    implementation("com.github.haifengl:smile-core:4.3.0")
    implementation("com.github.samtools:htsjdk:4.2.0")
    implementation("commons-cli:commons-cli:1.9.0")
    testImplementation(platform("org.junit:junit-bom:5.10.0"))
    testImplementation("org.junit.jupiter:junit-jupiter")
    testRuntimeOnly("org.junit.platform:junit-platform-launcher")
}

application {
    mainClass = "analysis.GenomeAnonymizer"
}

tasks.jar {
    manifest.attributes["Main-Class"] = "analysis.GenomeAnonymizer"
    val dependencies = configurations
        .runtimeClasspath
        .get()
        .map(::zipTree) // OR .map { zipTree(it) }
    from(dependencies)
    duplicatesStrategy = DuplicatesStrategy.EXCLUDE
}

// Unit tests only. Runs as part of ./gradlew build (skip with ./gradlew build -x test)
tasks.test {
    useJUnitPlatform {
        excludeTags("simulated")
    }
    // AnonymizedReadAlignmentProvider preallocates large collections
    maxHeapSize = "2g"
}

// Validation of the built JAR on the simulated tumor-normal datasets of src/test/resources/simulated
//   ./gradlew simulatedTest                    run the validation and compare against the committed baseline
//   ./gradlew simulatedTest -PupdateBaseline   accept the current results as the new baseline
//   ./gradlew simulatedTest -PstrictBaseline   fail on any difference against the baseline
val simulatedTest = tasks.register<Test>("simulatedTest") {
    description = "Runs the built JAR on the simulated tumor-normal datasets and validates the output"
    group = "verification"
    testClassesDirs = sourceSets.test.get().output.classesDirs
    classpath = sourceSets.test.get().runtimeClasspath
    useJUnitPlatform {
        includeTags("simulated")
    }
    dependsOn(tasks.jar)
    shouldRunAfter(tasks.test)
    systemProperty("simulated.jar", tasks.jar.get().archiveFile.get().asFile.absolutePath)
    systemProperty("simulated.dataDir", file("src/test/resources/simulated").absolutePath)
    systemProperty("simulated.outputDir", layout.buildDirectory.dir("simulated-test").get().asFile.absolutePath)
    systemProperty("simulated.reportDir", layout.buildDirectory.dir("reports/simulated").get().asFile.absolutePath)
    systemProperty("simulated.threads", (findProperty("simulatedThreads") ?: "8").toString())
    systemProperty("simulated.updateBaseline", project.hasProperty("updateBaseline").toString())
    systemProperty("simulated.strictBaseline", project.hasProperty("strictBaseline").toString())
    // The validation depends on the JAR and on files outside the task inputs, so it always runs when requested
    outputs.upToDateWhen { false }
    testLogging {
        events("passed", "failed", "skipped")
        showStandardStreams = true
        exceptionFormat = org.gradle.api.tasks.testing.logging.TestExceptionFormat.FULL
    }
}

// Build + unit tests + simulated validation
tasks.register("fullCheck") {
    description = "Builds the JAR and runs the unit tests and the simulated validation"
    group = "verification"
    dependsOn(tasks.build, simulatedTest)
}