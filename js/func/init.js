//global variables
var scene, camera, renderer, stats;

init();
animate();
function init() {
    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);

    camera.position.z = 5;

    renderer = new THREE.WebGLRenderer();
    renderer.setSize(window.innerWidth, window.innerHeight);
    document.body.appendChild(renderer.domElement);

    //stats fps
    stats = new Stats();
    var container = document.getElementById("fpsMeter");
    container.appendChild(stats.dom);
}
function animate() {
    requestAnimationFrame(animate);
    stats.update();
    renderer.render(scene, camera);
}