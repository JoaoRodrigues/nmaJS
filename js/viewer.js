// Simple Viewer
var glviewer = null;

export function displayPDB(moldata) {

    // Make gldiv and gradient legend visible and hide progress box
    let progressbox = document.getElementById("progress-text");
    progressbox.className = "hidden";

    let container = document.getElementById("viewer-container");
    container.className = "visible";

    // Now do actual stuff!
    glviewer = $3Dmol.createViewer("gldiv", {});

    // // resi-based gradient color scheme
    // let sortedResi = moldata.concat().sort(
    //     (a, b) => a.resi - b.resi
    // )

    // const colorByResi = {
    //     prop: 'resi',
    //     gradient: 'sinebow',
    //     min: sortedResi[0].resi,
    //     max: sortedResi[sortedResi.length - 1].resi
    // }

    let bfactors = moldata.map(bb => bb.b)
    const colorByFluct = {
        prop: 'b',
        gradient: 'rwb',
        min: Math.max(...bfactors),  // to reverse gradient
        max: Math.min(...bfactors)
    }

    let m = glviewer.addModel();
    m.addAtoms(moldata);

    // Assign callbacks to atoms
    const atoms = m.selectedAtoms({})
    for (let i in atoms) {
        let atom = atoms[i];
        atom.clickable = true;
        atom.callback = atomcallback;
    }

    // Add styles
    glviewer.setStyle({
        atom: ["CA", "P"],
    },
    {
        cartoon: {
            style: 'trace',
            ribbon: true,
            thickness: 1.0,
            colorscheme: colorByFluct,
        },
    }, )

    // glviewer.addStyle({
    //     atom: ["CA", "C4*", "COM"]
    // },
    // {
    //     stick: {
    //         radius: 0.5,
    //         colorscheme: colorByFluct
    //     },
    // });

//     glviewer.addStyle({
//         atom: ["COM"]
//     },
//     {
//         sphere: {
//             scale: 0.5,
//             colorscheme: colorByFluct
//         }
//     });

    glviewer.zoomTo();
    glviewer.render();
}


var atomcallback = function (atom, viewer) {
    if (atom.clickLabel === undefined ||
        !atom.clickLabel instanceof $3Dmol.Label) {
        atom.clickLabel = viewer.addLabel(atom.resn + atom.resi + '.' + atom.atom, {
            fontSize: 12,
            position: {
                x: atom.x,
                y: atom.y,
                z: atom.z
            },
            backgroundColor: "black"
        });
        atom.clicked = true;
    }

    //toggle label style
    else {

        if (atom.clicked) {
            var newstyle = atom.clickLabel.getStyle();
            newstyle.backgroundColor = 0x66ccff;

            viewer.setLabelStyle(atom.clickLabel, newstyle);
            atom.clicked = !atom.clicked;
        } else {
            viewer.removeLabel(atom.clickLabel);
            delete atom.clickLabel;
            atom.clicked = false;
        }

    }
};
