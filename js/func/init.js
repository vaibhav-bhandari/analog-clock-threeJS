//global variables
var scene, camera, renderer, stats, fontLoader;
var time, hrHand, minHand, secHand,
    handColor1 = 0x242424,
    handColor2 = 0x2762f3,
    mainColor = 0x242424,
    markingColor = 0x171717;

init();
animate();
function init() {
    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
    fontLoader = new THREE.FontLoader()

    time = new Time();
    time.update();

    var ambientLight = new THREE.AmbientLight("#ffffff");
    ambientLight.name = "Ambient Light";
    scene.add(ambientLight);

    var spotlightF = new THREE.SpotLight(0xffffff);
    spotlightF.name = "Front Spotlight";
    spotlightF.position.set(0, 20, 20);
    spotlightF.castShadow = true;
    spotlightF.shadow.mapSize = new THREE.Vector2(1024, 1024);
    spotlightF.shadow.camera.far = 40;
    spotlightF.shadow.camera.near = 0.5;
    scene.add(spotlightF);

    var spotlightB = new THREE.SpotLight(0xffffff);
    spotlightB.name = "Back Spotlight";
    spotlightB.position.set(0, 20, -20);
    spotlightB.castShadow = true;
    spotlightB.shadow.mapSize = new THREE.Vector2(1024, 1024);
    spotlightB.shadow.camera.far = 40;
    spotlightB.shadow.camera.near = 0.5;
    scene.add(spotlightB);

    var clock = new THREE.Object3D();
    clock.name = "Clock";
    scene.add(clock);

    var clockRimGeo = new THREE.TorusGeometry(10.5, 1, 16, 32),
        clockRimMat = new THREE.MeshStandardMaterial({
            color: mainColor
        }),
        clockRim = new THREE.Mesh(clockRimGeo, clockRimMat);
    clockRim.name = "Rim";
    clockRim.castShadow = true;
    clockRim.receiveShadow = true;
    clock.add(clockRim)

    var clockFaceGeo = new THREE.CylinderGeometry(10, 10, 0.5, 32, 32, false),
        clockFaceMat = new THREE.MeshStandardMaterial({
            color: 0xffffff
        }),
        clockFace = new THREE.Mesh(clockFaceGeo, clockFaceMat);
    clockFace.name = "Face";
    clockFace.position.set(0, 0, 0.25);
    clockFace.rotation.x = 90 * Math.PI / 180;
    clockFace.castShadow = true;
    clockFace.receiveShadow = true;
    clock.add(clockFace);

    fontLoader.load("https://threejs.org/examples/fonts/helvetiker_regular.typeface.json", (font) => {
        let fromCenter = 6.75,
            angleInc = 30;

        // 1. Digits
        for (let i = 1; i <= 12; ++i) {
            var digitGeo = new THREE.TextGeometry(`${i}`, {
                font: font,
                size: 2,
                height: 0,
                curveSegments: 16,
                bevelEnabled: false
            }),
                digitMat = new THREE.MeshStandardMaterial({
                    color: markingColor
                }),
                digit = new THREE.Mesh(digitGeo, digitMat);

            digit.name = `Digit ${i}`;
            digit.position.set(
                fromCenter * Math.sin((180 + -angleInc * i) * Math.PI / 180),
                0.26,
                fromCenter * Math.cos((180 + -angleInc * i) * Math.PI / 180)
            );
            digit.receiveShadow = true;
            digit.rotation.x = -90 * Math.PI / 180;

            digitGeo.computeBoundingBox();
            digitGeo.translate(-digitGeo.boundingBox.max.x / 2, -digitGeo.boundingBox.max.y / 2, 0);

            clockFace.add(digit);
        }
        var madeInIndiaLabel = "MADE IN INDIA",
            madeInIndiaGeo = new THREE.TextGeometry(madeInIndiaLabel, {
                font: font,
                size: 0.25,
                height: 0.1,
                curveSegments: 16,
                bevelEnabled: false
            }),
            madeInIndiaMat = new THREE.MeshStandardMaterial({
                color: mainColor
            }),
            madeInIndia = new THREE.Mesh(madeInIndiaGeo, madeInIndiaMat);
        madeInIndia.name = madeInIndiaLabel;
        madeInIndia.position.set(0.3, -2.7, -2.15);
        madeInIndia.rotation.x = 180 * Math.PI / 180;
        madeInIndia.rotation.z = 180 * Math.PI / 180;
        madeInIndia.receiveShadow = true;
        madeInIndiaGeo.computeBoundingBox();
        clock.add(madeInIndia);
    });

    var thickTickGeo = new THREE.BoxGeometry(0.3, 0.1, 1),
        thickTickMat = new THREE.MeshStandardMaterial({
            color: markingColor
        }),
        thickTick = new THREE.Mesh(thickTickGeo, thickTickMat);
    thickTick.receiveShadow = true;

    var thinTickGeo = new THREE.BoxGeometry(0.15, 0.1, 0.75),
        thinTickMat = new THREE.MeshStandardMaterial({
            color: markingColor
        }),
        thinTick = new THREE.Mesh(thinTickGeo, thinTickMat);
    thinTick.receiveShadow = true

    let angleInc = 6;
    for (let i = 0; i < 60; ++i) {
        let divBy5 = !(i % 5),
            tick = divBy5 ? thickTick.clone() : thinTick.clone(),
            fromCenter = divBy5 ? 8.875 : 9;

        tick.name = `Tick ${i}`;
        tick.position.set(
            fromCenter * Math.sin((-angleInc * i) * Math.PI / 180),
            0.31,
            fromCenter * Math.cos((-angleInc * i) * Math.PI / 180)
        );
        tick.rotation.y = (-angleInc * i) * Math.PI / 180;
        clockFace.add(tick);
    }

    // E. Hands
    hrHand = new Hand({
        group: clock,
        w: 0.5,
        h: 5,
        d: 0.25,
        botExt: 2,
        z: 0.85,
        rot: time.getHrAngle(),
        color: handColor1,
        name: "Hour Hand"
    });
    minHand = new Hand({
        group: clock,
        w: 0.5,
        h: 8.5,
        d: 0.25,
        botExt: 2,
        z: 1.1,
        rot: time.getMinAngle(),
        moveDivs: 60 * 60,
        color: handColor1,
        name: "Minute Hand"
    });
    secHand = new Hand({
        group: clock,
        w: 0.25,
        h: 8.5,
        d: 0.25,
        botExt: 2,
        z: 1.35,
        rot: time.getSecAngle(),
        moveDivs: 60,
        color: handColor2,
        showHub: true,
        name: "Second Hand"
    });

    // F. Dome
    var domeGeo = new THREE.SphereGeometry(9.625, 32, 8, 0, Math.PI * 2, 0, Math.PI / 2),
        domeMat = new THREE.MeshStandardMaterial({
            color: 0xffffff,
            opacity: 0.3,
            transparent: true
        }),
        dome = new THREE.Mesh(domeGeo, domeMat);
    dome.position.y = 0.25;
    dome.scale.y = 0.25;
    clockFace.add(dome);

    // G. Back
    // 1. Main
    var clockBackGeo = new THREE.CylinderGeometry(10, 10, 0.5, 32, 32, false),
        clockBackMat = new THREE.MeshStandardMaterial({
            color: mainColor
        }),
        clockBack = new THREE.Mesh(clockBackGeo, clockBackMat);
    clockBack.name = "Back";
    clockBack.position.set(0, 0, -0.25);
    clockBack.rotation.x = 90 * Math.PI / 180;
    clockBack.castShadow = true;
    clockBack.receiveShadow = true;
    clock.add(clockBack);

    // 2. Back Box
    var backBoxBase = new THREE.Mesh(new THREE.BoxGeometry(6, 1.5, 7)),
        nailHoleA = new THREE.Mesh(new THREE.CylinderGeometry(0.5, 0.5, 1, 16, 16, false)),
        nailHoleB = new THREE.Mesh(new THREE.CylinderGeometry(0.3, 0.25, 1, 16, 16, false));

    nailHoleA.position.add(new THREE.Vector3(0, -0.5, -2));
    nailHoleB.position.add(new THREE.Vector3(0, -0.5, -2.5));

    backBoxBase.updateMatrix();
    nailHoleA.updateMatrix();
    nailHoleB.updateMatrix();

    backBoxBase.material = new THREE.MeshStandardMaterial({
        color: mainColor
    });

    //cut holes for the nail, mix it altogether
    // var backBoxBase_BSP = CSG.fromMesh(backBoxBase),
    //     nailHoleA_BSP = CSG.fromMesh(nailHoleA),
    //     nailHoleB_BSP = CSG.fromMesh(nailHoleB),
    //     backBox_BSP = backBoxBase_BSP.subtract(nailHoleA_BSP).subtract(nailHoleB_BSP),
    //     backBox = CSG.toMesh(backBox_BSP, backBoxBase.matrix, backBoxBase.material);

    
    // backBox.name = "Back Box";
    // backBox.position.set(0, -1, 0);
    // backBox.castShadow = true;
    // backBox.receiveShadow = true;
    // clockBack.add(backBox);

    // 3. Battery Door
    var batteryDoorGeo = new THREE.BoxGeometry(5, 0.2, 3.5),
        batteryDoorMat = new THREE.MeshStandardMaterial({
            color: mainColor
        }),
        batteryDoor = new THREE.Mesh(batteryDoorGeo, batteryDoorMat);
    batteryDoor.name = "Battery Door";
    batteryDoor.position.set(0, -0.85, 1.25);
    batteryDoor.castShadow = true;
    batteryDoor.receiveShadow = true;
    //backBox.add(batteryDoor);

    // IV. Camera
    camera.position.set(0, 0, 45);
    camera.lookAt(scene.position);

    renderer = new THREE.WebGLRenderer();
    renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.setClearColor(new THREE.Color(0x87a9f9));
    renderer.shadowMap.enabled = true;
    document.body.appendChild(renderer.domElement);

    //stats fps
    stats = new Stats();
    var container = document.getElementById("fpsMeter");
    container.appendChild(stats.dom);
    let camControls = new THREE.OrbitControls(camera, renderer.domElement);
    
}
function animate() {
    requestAnimationFrame(animate);
    time.update();
    hrHand.rotate(time.getHrAngle());
    minHand.rotate(time.getMinAngle(), 0.1 / 6);
    secHand.rotate(time.getSecAngle(), 1);
    stats.update();
    renderer.render(scene, camera);
}

adjustScreen = () => {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight)
};
window.addEventListener("resize",adjustScreen);