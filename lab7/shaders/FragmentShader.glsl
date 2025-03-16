#version 330 core

in vec3 color;
out vec4 out_color;

uniform float altitude_player;
uniform float max_altitude;

uniform bool is_terrain;

void main()
{
    if (is_terrain) {
        vec3 color1 = vec3(0.55, 0.27, 0.07);  // Maro

        vec3 color2 = vec3(0.0, 1.0, 0.0);  // Verde

        // Interpolare între cele două culori în funcție de altitudine
        float t = clamp(altitude_player / max_altitude, 0.0, 1.0);
        vec3 interpolatedColor = mix(color1, color2, t);

       // out_color = vec4(interpolatedColor, 1.0);
       out_color = vec4(color, 1.0);
    } else {
		//culoarea initiala
        out_color = vec4(color, 1.0);
    }
}