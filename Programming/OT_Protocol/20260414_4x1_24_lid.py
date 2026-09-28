from opentrons import protocol_api, types

# metadata
metadata = {
    "author": "Mahdi Rastegar - Modified for Release Experiment with Lid Handling",
    "apiLevel": "2.18"
}

# Protocol for Release Experiment (single apparatus) with automated lid handling
# 1- manually loadings at lower row at HTDA block
# 2- OT wait 5 minutes before starting from the moment we loaded sample manually
# 3- For every time point:
#      a. OT removes the lid from the wellplate (slot 3) and parks it on the lid holder (slot 6)
#      b. OT will take sample from upper row and transfer to Microplate 250 ul (all 4 cells)
#      c. OT removes the lid from the holder and places it back on the wellplate
#      d. OT replace media with 250 ul from reservoir to upper row
# 4- OT will repeat step 3 for 24 times with 12 minutes interval (using rows A-D for first 12, E-H for next 12)

# protocol run function
def run(protocol: protocol_api.ProtocolContext):

    # labware
    release_experiment_1 = protocol.load_labware("sdl5_4_channel_apparatus_12ul", location="4")

    tiprack_1 = protocol.load_labware("opentrons_96_tiprack_300ul", location="2")
    # tiprack_2 = protocol.load_labware("opentrons_96_tiprack_300ul", location="5")

    release_sample = protocol.load_labware("greiner_96_wellplate_382ul", location="3")
    release_media = protocol.load_labware("calab_8_tuberack_20000ul", location="7")

    # Lid holder for the sample wellplate's lid
    lid_holder = protocol.load_labware("greiner_96_wellplate_382ul", location="6")  # Assuming lid is stored/parked in a 96-well plate format for simplicity

    # Denote an empty slot as a 1000 uL tiprack - used as a "fake" tip so the lid tool can grip/release the lid
    tiprack_fake = protocol.load_labware("opentrons_lid_tiprack_1000ul", location="11")

    # pipettes
    p300 = protocol.load_instrument("p300_single", mount="left", tip_racks=[tiprack_1])
    # p300 = protocol.load_instrument("p300_single", mount="left", tip_racks=[tiprack_1, tiprack_2])

    # Right pipette: mechanical lid handling via pick_up_tip/drop_tip at calibrated points.
    lid_tool = protocol.load_instrument("p1000_single", mount="right")

    # experiment setup
    sample_volume = 250  # volume to transfer from upper row to microplate
    media_volume = 250   # volume to replace in upper row from reservoir
    p300.flow_rate.aspirate = 50
    p300.flow_rate.dispense = 50

    # Lid holder and plate are treated like a 96-well grid; D6 + XY offset hits physical lid center.
    lid_center_well = "D6"
    lid_center_xy_offset = (0, 0)
    lid_center_z_offset = 20
    # Grip (pick up) depth is confirmed working. Release depth is separately tunable:
    # a friction/collet fit often needs to be pushed measurably deeper to let go than
    # to seat. If release still doesn't drop the lid, jog deeper in ~2-5mm increments.
    lid_grip_depth = 18.8
    lid_release_depth = 18.8

    def take_sample_from_upper_row(sample_volume, time_point):
        """Take sample from upper row (row A) and transfer to microplate"""
        p300.well_bottom_clearance.dispense = 5
        p300.well_bottom_clearance.aspirate = 2

        # Each apparatus cell goes to its corresponding wellplate rows, different columns for each time point
        # Time points 0-11: use rows A,B,C,D; Time points 12-23: use rows E,F,G,H
        if time_point < 12:
            microplate_rows = ['A', 'B', 'C', 'D']
            column_number = time_point + 1  # Columns 1-12 for first set
        else:
            microplate_rows = ['E', 'F', 'G', 'H']
            column_number = (time_point - 12) + 1  # Columns 1-12 for second set

        for n in range(4):  # 4 wells A1-A4 from apparatus
            p300.pick_up_tip()
            from_well = 'A' + str(n+1)  # Upper row wells A1, A2, A3, A4 from apparatus
            target_row = microplate_rows[n]  # A1→A, A2→B, A3→C, A4→D
            to_well = target_row + str(column_number)  # A1-A12,B1-B12,C1-C12,D1-D12, then E1-E12,F1-F12,G1-G12,H1-H12
            p300.transfer(
                sample_volume,
                release_experiment_1[from_well],
                release_sample[to_well],
                blow_out=True,
                blowout_location='destination well',
                new_tip='never'
            )
            p300.return_tip()

    def replace_media_in_upper_row(media_volume):
        """Replace media with 250 ul from reservoir to upper row"""
        p300.well_bottom_clearance.dispense = 3
        p300.well_bottom_clearance.aspirate = 3

        for n in range(4):  # 4 wells A1-A4
            p300.pick_up_tip()
            from_well = 'A' + str(n+1)  # Reservoir wells A1, A2, A3, A4
            to_well = 'A' + str(n+1)    # Upper row wells A1, A2, A3, A4
            p300.transfer(
                media_volume,
                release_media[from_well],
                release_experiment_1[to_well],
                blow_out=True,
                blowout_location='destination well',
                new_tip='never'
            )
            p300.return_tip()

    def lid_top_with_xy_offset(target_well, z_height, xy_offset):
        # Build a centered lid contact location from nominal well top + calibrated XY shift.
        return target_well.top(z=z_height).move(types.Point(x=xy_offset[0], y=xy_offset[1], z=0))

    def custom_lid_approach(instr, target_well, z_offset=lid_center_z_offset, approach_z=8, lid_depth=lid_grip_depth, xy_offset=(5.0, -5.0), dwell=0.5):
        # Mimics the old custom sequence: approach above lid, descend to grip/release depth, retract.
        instr.move_to(lid_top_with_xy_offset(target_well, z_offset + approach_z, xy_offset))
        protocol.delay(seconds=dwell)
        instr.move_to(lid_top_with_xy_offset(target_well, z_offset - lid_depth, xy_offset))
        protocol.delay(seconds=dwell)
        instr.move_to(lid_top_with_xy_offset(target_well, z_offset, xy_offset))
        protocol.delay(seconds=dwell)

    def remove_lid_from_wellplate():
        """Move the lid from the sample wellplate (slot 3) onto the lid holder (slot 6)"""
        plate_center = release_sample[lid_center_well]
        holder_well = lid_holder[lid_center_well]

        # 1. Pick up the well plate lid with the motion to physically pick up the lid (microplate)
        custom_lid_approach(lid_tool, plate_center, z_offset=lid_center_z_offset, lid_depth=lid_grip_depth, xy_offset=lid_center_xy_offset)
        lid_tool.move_to(lid_top_with_xy_offset(plate_center, lid_center_z_offset - 60, lid_center_xy_offset))

        # 2. Pick up an "imaginary" pipette tip to tell the OT2 it has a tip (which is actually our lid)
        lid_tool.pick_up_tip(location=tiprack_fake["D6"])  # pick up from the middle to minimize chances for collision

        # 3. Complete the drop off of the lid onto the holder plate
        lid_tool.move_to(lid_top_with_xy_offset(holder_well, lid_center_z_offset - 60, lid_center_xy_offset))
        custom_lid_approach(lid_tool, holder_well, z_offset=lid_center_z_offset, lid_depth=lid_release_depth, xy_offset=lid_center_xy_offset)
        lid_tool.move_to(lid_top_with_xy_offset(holder_well, lid_center_z_offset, lid_center_xy_offset))

        # 4. Drop the lid onto the holder well, if needed, change drop height if too high.
        lid_tool.drop_tip(location=holder_well.top(-30))

    def place_lid_on_wellplate():
        """Move the lid from the lid holder (slot 6) back onto the sample wellplate (slot 3)"""
        plate_center = release_sample[lid_center_well]
        holder_well = lid_holder[lid_center_well]

        # 1. Pick up the well plate lid with the motion to physically pick up the lid (holder)
        custom_lid_approach(lid_tool, holder_well, z_offset=lid_center_z_offset, lid_depth=lid_grip_depth, xy_offset=lid_center_xy_offset)
        lid_tool.move_to(lid_top_with_xy_offset(holder_well, lid_center_z_offset - 40, lid_center_xy_offset))

        # 2. Pick up an "imaginary" pipette tip to tell the OT2 it has a tip (which is actually our lid)
        lid_tool.pick_up_tip(location=tiprack_fake["D6"])  # pick up from the middle to minimize chances for collision

        # 3. Complete the drop off of the lid onto the microplate
        lid_tool.move_to(lid_top_with_xy_offset(plate_center, lid_center_z_offset, lid_center_xy_offset))
        custom_lid_approach(lid_tool, plate_center, z_offset=lid_center_z_offset, lid_depth=lid_release_depth, xy_offset=lid_center_xy_offset)
        lid_tool.move_to(lid_top_with_xy_offset(plate_center, lid_center_z_offset - 20, lid_center_xy_offset))

        # 4. Drop the lid onto the plate, if needed, change drop height if too high.
        lid_tool.drop_tip(location=plate_center.top(-30))

    # Run the protocol
    protocol.comment("Step 1: Manual loading completed at lower row of HTDA block")
    protocol.comment("Step 2: Waiting 5 minutes before starting automated protocol")

    # Wait 5 minutes for manual loading
    protocol.delay(minutes=5)

    # Repeat sampling and media replacement 24 times with 12-minute intervals
    for time_point in range(24):
        # Reset tip racks after 12 time points to reuse tips for the next 12 time points
        #--------------------------------------#
        if time_point == 12:
            p300.reset_tipracks()
            protocol.comment("Tip rack reset: reusing tiprack in location 2 for time points 13-24")
        #--------------------------------------#
        protocol.comment(f"Starting time point {time_point + 1} of 24")

        # Step 3a: Remove lid from wellplate and park it on the lid holder
        protocol.comment("Step 3a: Removing lid from wellplate, parking on lid holder (slot 6)")
        remove_lid_from_wellplate()

        # Step 3b: Take sample from upper row and transfer to microplate
        protocol.comment(f"Step 3b: Taking {sample_volume} µL sample from upper row to microplate")
        take_sample_from_upper_row(sample_volume, time_point)

        # Step 3c: Place lid back on wellplate
        protocol.comment("Step 3c: Placing lid back on wellplate")
        place_lid_on_wellplate()

        # Step 3d: Replace media with 250 µL from reservoir to upper row
        protocol.comment(f"Step 3d: Replacing with {media_volume} µL fresh media from reservoir")
        replace_media_in_upper_row(media_volume)

        # Wait 12 minutes before next cycle (except after the last cycle)
        if time_point < 23:  # Don't wait after the 24th and final time point
            protocol.comment("Waiting 12 minutes before next sampling time point")
            protocol.delay(minutes=12)

    protocol.comment("Release experiment protocol completed!")
