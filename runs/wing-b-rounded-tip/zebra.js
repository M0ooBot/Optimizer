function createZebraRenderer(canvas, data, cameraPresets) {
    const gl = canvas.getContext("webgl2", {antialias: true, alpha: false});
    if (!gl) throw new Error("WebGL 2 is unavailable");

    const vertex = `#version 300 es
        in vec3 aPosition;
        in vec3 aNormal;
        uniform vec3 uCenter, uAcross, uUp, uView;
        uniform vec3 uScale;
        uniform float uSide;
        out vec3 normal;
        out float span;
        void main() {
            vec3 mirror = vec3(uSide, 1.0, 1.0);
            vec3 position = aPosition * mirror - uCenter;
            normal = aNormal * mirror;
            span = aPosition.x;
            gl_Position = vec4(dot(position,uAcross)*uScale.x, dot(position,uUp)*uScale.y,
                               -dot(position,uView)*uScale.z, 1.0);
        }`;
    const fragment = `#version 300 es
        precision highp float;
        in vec3 normal;
        in float span;
        uniform vec3 uView, uUp, uAcross, uBandNormal, uSurface;
        uniform float uStart, uAlpha, uSide;
        uniform bool uZebra;
        out vec4 stripeColor;
        void main() {
            if (span < uStart) discard;
            vec3 n = normalize(normal);
            if (!uZebra) {
                vec3 lightDir = normalize(0.34*uAcross + 0.79*uUp + 0.55*uView);
                float diffuse = max(dot(n, lightDir), 0.0);
                float sheen = pow(max(dot(n, normalize(lightDir + uView)), 0.0), 27.0) * 0.18;
                stripeColor = vec4(uSurface * (0.57 + 0.43 * diffuse) + sheen, uAlpha);
                return;
            }
            // Reflect the strip-light rig across the root with the left wing.
            n.x *= uSide;
            vec3 incident = vec3(-uView.x * uSide, -uView.y, -uView.z);
            vec3 reflection = reflect(incident, n);
            float phase = 22.0 * dot(reflection, uBandNormal);
            float footprint = fwidth(phase);
            float edge = clamp(footprint * 0.47, 0.025, 0.95);
            float band = smoothstep(-edge, edge, sin(phase));
            band = mix(band, 0.5, smoothstep(0.85, 3.1, footprint));
            vec3 dark = vec3(0.045, 0.07, 0.075);
            vec3 light = vec3(0.925, 0.942, 0.937);
            stripeColor = vec4(mix(dark, light, band), 1.0);
        }`;

    function compile(type, source) {
        const shader = gl.createShader(type);
        gl.shaderSource(shader, source);
        gl.compileShader(shader);
        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(shader));
        return shader;
    }
    const program = gl.createProgram();
    gl.attachShader(program, compile(gl.VERTEX_SHADER, vertex));
    gl.attachShader(program, compile(gl.FRAGMENT_SHADER, fragment));
    gl.linkProgram(program);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(program));
    gl.useProgram(program);
    gl.enable(gl.DEPTH_TEST);

    const uniforms = Object.fromEntries(
        ["Center", "Across", "Up", "View", "Scale", "Side", "BandNormal", "Start",
         "Surface", "Alpha", "Zebra"]
            .map(name => [name, gl.getUniformLocation(program, "u" + name)])
    );
    function typedArray(encoded, ArrayType) {
        const encodedBytes = atob(encoded);
        const bytes = Uint8Array.from(encodedBytes, char => char.charCodeAt(0));
        return new ArrayType(bytes.buffer);
    }
    function upload(target, array) {
        const buffer = gl.createBuffer();
        gl.bindBuffer(target, buffer);
        gl.bufferData(target, array, gl.STATIC_DRAW);
        return buffer;
    }
    function attribute(name, encoded) {
        upload(gl.ARRAY_BUFFER, typedArray(encoded, Float32Array));
        const index = gl.getAttribLocation(program, name);
        gl.vertexAttribPointer(index, 3, gl.FLOAT, false, 0, 0);
        gl.enableVertexAttribArray(index);
    }
    function indices(encoded) {
        const array = typedArray(encoded, Uint32Array);
        return {buffer: upload(gl.ELEMENT_ARRAY_BUFFER, array), count: array.length};
    }
    const meshes = Object.fromEntries(Object.entries(data).map(([name, mesh]) => {
        const vao = gl.createVertexArray();
        gl.bindVertexArray(vao);
        attribute("aPosition", mesh.positions);
        attribute("aNormal", mesh.normals);
        return [name, {vao, solid: indices(mesh.triangles), skin: indices(mesh.skin)}];
    }));
    gl.bindVertexArray(null);

    const bounds = {
        "Full wing": {center: [0, 100, 31], size: [1200, 200, 72]},
        "Full half-wing": {center: [300, 100, 31], size: [600, 200, 72]},
        "Tip detail": {center: [563, 110, 58], size: [76, 97, 22]},
    };
    const byCoordinate = ({x, y, z}) => [x, y, z];
    const norm = vector => vector.map(value => value / Math.hypot(...vector));
    function cross(a, b) {
        return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]];
    }
    function basis() {
        const preset = cameraPresets[state.scope === "Full wing" && state.camera === "Angled"
            ? "Full angled" : state.camera];
        if (!orbit) {
            const view = norm(byCoordinate(preset.eye));
            const across = norm(cross(byCoordinate(preset.up), view));
            return [view, across, cross(view, across)];
        }
        const [azimuth, elevation] = orbit;
        const view = [Math.cos(elevation)*Math.cos(azimuth),
                      Math.cos(elevation)*Math.sin(azimuth), Math.sin(elevation)];
        const across = norm(cross([0, 0, 1], view));
        return [view, across, cross(view, across)];
    }
    let state;
    let orbit = null;
    let zoom = 1;
    let pointer;

    function paint() {
        if (!state || canvas.hidden || !canvas.clientWidth) return;
        const density = Math.min(devicePixelRatio || 1, 2);
        const width = Math.round(canvas.clientWidth * density);
        const height = Math.round(canvas.clientHeight * density);
        if (canvas.width !== width) canvas.width = width;
        if (canvas.height !== height) canvas.height = height;
        const colors = ["Original", "Rounded tip"];
        const compare = state.model === "Overlay";
        const zebra = state.shading === "Zebra";
        const sideBySide = compare && zebra;
        const models = compare ? colors : [state.model];
        const [view, across, up] = basis();
        const region = bounds[state.scope];
        const full = state.scope === "Full wing";
        const extent = axis => axis.reduce((total, n, i) => total + Math.abs(n)*region.size[i], 0);
        const angle = state.rotation * Math.PI / 180;
        gl.clearColor(.997, .998, 1, 1);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        gl.uniform3fv(uniforms.View, view);
        gl.uniform3fv(uniforms.Across, across);
        gl.uniform3fv(uniforms.Up, up);
        gl.uniform3fv(uniforms.Center, region.center);
        gl.uniform3f(uniforms.BandNormal, Math.sin(angle), 0, Math.cos(angle));
        gl.uniform1f(uniforms.Start, state.scope === "Tip detail" ? 525 : -1);
        gl.uniform1i(uniforms.Zebra, zebra);
        models.forEach((model, index) => {
            const left = sideBySide ? Math.floor(index * width / 2) : 0;
            const right = sideBySide ? Math.floor((index + 1) * width / 2) : width;
            const aspect = (right - left) / height;
            const span = Math.max(extent(across), extent(up) * aspect) / (full ? .88 : .86);
            gl.viewport(left, 0, right - left, height);
            gl.uniform3f(uniforms.Scale, 2*zoom/span, 2*zoom*aspect/span, 1/region.size[0]);
            const mesh = meshes[model];
            gl.bindVertexArray(mesh.vao);
            const draw = full ? mesh.skin : mesh.solid;
            gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, draw.buffer);
            gl.uniform3fv(uniforms.Surface, model === "Original" ? [.38, .48, .54] : [.13, .64, .61]);
            gl.uniform1f(uniforms.Alpha, compare && !zebra ? .46 : 1);
            function drawModel() {
                gl.uniform1f(uniforms.Side, 1);
                gl.drawElements(gl.TRIANGLES, draw.count, gl.UNSIGNED_INT, 0);
                if (full) {
                    gl.uniform1f(uniforms.Side, -1);
                    gl.drawElements(gl.TRIANGLES, draw.count, gl.UNSIGNED_INT, 0);
                }
            }
            if (compare && !zebra) {
                // Blend each model only where its own visible surface passed the depth test.
                gl.clear(gl.DEPTH_BUFFER_BIT);
                gl.colorMask(false, false, false, false);
                drawModel();
                gl.colorMask(true, true, true, true);
                gl.depthFunc(gl.EQUAL);
                gl.enable(gl.BLEND);
                gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
                drawModel();
                gl.disable(gl.BLEND);
                gl.depthFunc(gl.LESS);
            } else {
                drawModel();
            }
        });
        gl.bindVertexArray(null);
    }

    canvas.addEventListener("pointerdown", event => {
        if (!state) return;
        const [view] = basis();
        if (!orbit) orbit = [Math.atan2(view[1], view[0]), Math.asin(view[2])];
        pointer = [event.clientX, event.clientY];
        canvas.setPointerCapture(event.pointerId);
    });
    canvas.addEventListener("pointermove", event => {
        if (!pointer) return;
        orbit[0] -= (event.clientX - pointer[0]) * .008;
        orbit[1] = Math.max(-1.56, Math.min(1.56, orbit[1] + (event.clientY - pointer[1]) * .008));
        pointer = [event.clientX, event.clientY];
        paint();
    });
    canvas.addEventListener("lostpointercapture", () => { pointer = null; });
    canvas.addEventListener("wheel", event => {
        event.preventDefault();
        zoom = Math.max(.55, Math.min(9, zoom * Math.exp(-event.deltaY * .001)));
        paint();
    }, {passive: false});
    new ResizeObserver(paint).observe(canvas);
    function update(next) {
        if (!state || next.scope !== state.scope || next.camera !== state.camera || next.reset) {
            orbit = null;
            zoom = 1;
        }
        state = next;
        paint();
    }
    update.snapshot = () => {
        paint();
        return canvas.toDataURL("image/png");
    };
    return update;
}
