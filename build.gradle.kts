plugins {
    id("java")
    application
}

group = "org.example"
version = "1.0.0"

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

tasks.test {
    useJUnitPlatform()
}