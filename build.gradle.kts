plugins {
    id("java")
    application
}

group = "org.example"
version = "0.0.3"

repositories {
    mavenCentral()
}

dependencies {
    implementation(fileTree(mapOf("dir" to "lib", "include" to listOf("*.jar"))))
    // implementation("com.github.samtools:htsjdk:4.1.3")
    // implementation("commons-cli:commons-cli:1.9.0")
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