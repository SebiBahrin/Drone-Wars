#version 330 core

layout(location = 0) in vec3 v_position;
layout(location = 1) in vec3 v_normal;
layout(location = 2) in vec2 v_texture_coord;

uniform mat4 Model;
uniform mat4 View;
uniform mat4 Projection;

uniform bool is_terrain;
uniform float noise_frequency;
uniform float noise_amplitude;

uniform vec3 object_color;

out vec3 color;

float noise(vec2 st) {
    return fract(sin(dot(st.xy, vec2(12.9898, 78.233))) * 43758.5453123);
}

void main()
{
    vec3 position = v_position;

    // Aplicați noise doar pentru teren
    if (is_terrain) {
        // Calculăm noise-ul pe baza poziției (x, z) a vertecșilor
        float elevation = noise(position.xz * noise_frequency) * noise_amplitude;

        // Modificăm poziția pe axa Y folosind valoarea noise
        position.y += elevation;

        // Ajustăm culoarea în funcție de noise
        //color = mix(vec3(0.55, 0.27, 0.07), vec3(0.0, 1.0, 0.0), elevation);
        color = object_color;
    } else {
        // Pentru alte obiecte folosim culoarea transmisă din CPU
        color = object_color;
    }

    gl_Position = Projection * View * Model * vec4(position, 1.0);
}