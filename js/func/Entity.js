class Hand {
    constructor(args) {
        try {
            let group = args.group,
                w = args.w || 1,
                h = args.h || 10,
                d = args.d || 1,
                botExt = args.botExt || 0,
                z = args.z || 0,
                color = args.color || 0x000000,
                showHub = args.showHub === true ? true : false,
                name = args.name || "";

            // main
            this.rot = args.rot || 0;
            this.moveDivs = args.moveDivs || 1;

            // model
            this.geo = new THREE.BoxGeometry(w,h + botExt,d);
            this.mat = new THREE.MeshStandardMaterial({
                color: color
            });
            this.mesh = new THREE.Mesh(this.geo,this.mat);
            this.mesh.name = name;
            this.mesh.position.set(0,0,z);
            this.mesh.castShadow = true;
            this.mesh.receiveShadow = true;
            this.geo.translate(0,(h + botExt)/2 - botExt,0);
            group.add(this.mesh);

            if (showHub) {
                var hubGeo = new THREE.CylinderGeometry(w*2,w*2,d,16,16,false),
                    hubMat = this.mat,
                    hub = new THREE.Mesh(hubGeo,hubMat);
                hub.rotation.x = 90 * Math.PI/180;
                hub.castShadow = true;
                hub.receiveShadow = true;
                this.mesh.add(hub);
            }

        } catch {
            console.log("Object needs a parent group before adding to the scene.");
        }
    }
    rotate(deg,speed = 0) {
        try {
            // trick ensure proper transition going back to 0 tick from left
            if (deg == 0 && this.rot > 0) {
                this.rot = -360/this.moveDivs;
            }
            // normal transition
            if (this.rot < deg && speed > 0) {
                this.rot += speed;

            } else {
                this.rot = deg;
            }
            // rotate the shape
            this.mesh.rotation.z = -this.rot * Math.PI/180;

        } catch {
            console.log("You can’t rotate a nonexistent object.");
        }
    }
}

class Time {
    constructor() {
        this.stamp = null;
        this.hr = 0;
        this.min = 0;
        this.sec = 0;
    }
    update() {
        this.stamp = new Date();

        this.hr = this.stamp.getHours();
        if (this.hr > 12)
            this.hr -= 12;

        this.min = this.stamp.getMinutes();
        this.sec = this.stamp.getSeconds();
    }
    getHrAngle() {
        return (time.hr + time.min/60 + time.sec/(60*60)) / 12 * 360;
    }
    getMinAngle() {
        return (time.min + time.sec/60) / 60 * 360;
    }
    getSecAngle() {
        return time.sec/60 * 360;
    }
}